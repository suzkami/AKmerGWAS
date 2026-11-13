"""PLINK 1 ``.bed`` / ``.bim`` / ``.fam`` primitives.

Only *three* operations are needed by the pipeline:

* **encode** a presence/absence row directly into ``.bed`` bytes (the streaming
  multibed path, so no ``.tped`` ever touches disk),
* **concatenate** several chunk beds that share an identical sample order into
  one bed (used to build ``QK_bed`` and the per-trait ``siggeno`` bed without
  ``plink --merge-list``, which would silently reorder samples),
* **subset** a bed by variant id (used to pull significant k-mers out of a
  chunk).

Genotype encoding (PLINK 1, SNP-major):

    magic      = 6C 1B 01
    per SNP    = ceil(N/4) bytes
    code (2 bits) -> 0b00 hom for allele1, 0b01 missing, 0b10 het, 0b11 hom for allele2

Individual ``i`` occupies the two bits ``(2*(i%4))`` and ``(2*(i%4)+1)`` of byte
``i//4`` (the first sample is in the least-significant two bits). Because a bed
is a *pure concatenation* of independently encoded SNPs, the concatenate/subset
paths are bit-order agnostic -- they copy bytes verbatim.
"""

from __future__ import annotations

import os
from typing import Iterable, Iterator, Sequence

MAGIC = b"\x6c\x1b\x01"

# genotype codes
HOM_A1 = 0b00
MISSING = 0b01
HET = 0b10
HOM_A2 = 0b11


def bytes_per_variant(n_samples: int) -> int:
    return (n_samples + 3) // 4


def pack_variant(codes: Sequence[int]) -> bytes:
    """Pack per-sample genotype codes into ``ceil(n/4)`` bytes (LSB-first)."""
    buf = bytearray(bytes_per_variant(len(codes)))
    for i, code in enumerate(codes):
        buf[i >> 2] |= (code & 0x03) << ((i & 3) * 2)
    return bytes(buf)


def unpack_variant(buf: bytes, n_samples: int) -> list[int]:
    """Inverse of :func:`pack_variant` (used by the tests)."""
    return [(buf[i >> 2] >> ((i & 3) * 2)) & 0x03 for i in range(n_samples)]


def bim_line(chrom: str, vid: str, cm: str, pos: str, allele1: str, allele2: str) -> str:
    return "\t".join([chrom, vid, cm, pos, allele1, allele2])


# ---------------------------------------------------------------------------
# Reading headers
# ---------------------------------------------------------------------------

def read_fam(fam_path: str) -> list[str]:
    with open(fam_path, encoding="utf-8") as fh:
        return [ln.rstrip("\n") for ln in fh if ln.strip()]


def fam_sample_ids(fam_path: str) -> list[str]:
    """Second column (IID) of a ``.fam`` file, in file order."""
    ids = []
    for line in read_fam(fam_path):
        parts = line.split()
        ids.append(parts[1] if len(parts) > 1 else parts[0])
    return ids


def read_bim(bim_path: str) -> list[str]:
    with open(bim_path, encoding="utf-8") as fh:
        return [ln.rstrip("\n") for ln in fh if ln.strip()]


def iter_bim_ids(bim_path: str) -> Iterator[str]:
    with open(bim_path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.split()
            if parts:
                yield parts[1]


# ---------------------------------------------------------------------------
# Writing
# ---------------------------------------------------------------------------

class BedWriter:
    """Incrementally writes a ``.bed``/``.bim``/``.fam`` triple.

    Writes go to ``*.tmp`` siblings and are promoted with ``os.replace`` by
    :meth:`close`, so a crash never leaves a half-written chunk that later looks
    complete.
    """

    def __init__(self, prefix: str, fam_lines: Sequence[str]) -> None:
        self.prefix = prefix
        self.fam_lines = list(fam_lines)
        os.makedirs(os.path.dirname(os.path.abspath(prefix)) or ".", exist_ok=True)
        self._bed = open(prefix + ".bed.tmp", "wb")
        self._bed.write(MAGIC)
        self._bim = open(prefix + ".bim.tmp", "w", encoding="utf-8")
        self.n_variants = 0

    def add_variant(self, bim_fields: Sequence[str], packed: bytes) -> None:
        self._bim.write("\t".join(bim_fields) + "\n")
        self._bed.write(packed)
        self.n_variants += 1

    def close(self) -> None:
        self._bed.close()
        self._bim.close()
        with open(self.prefix + ".fam.tmp", "w", encoding="utf-8") as fh:
            for line in self.fam_lines:
                fh.write(line + "\n")
        for ext in ("bed", "bim", "fam"):
            os.replace(f"{self.prefix}.{ext}.tmp", f"{self.prefix}.{ext}")

    def abort(self) -> None:
        try:
            self._bed.close()
        except Exception:  # noqa: BLE001
            pass
        try:
            self._bim.close()
        except Exception:  # noqa: BLE001
            pass
        for ext in ("bed", "bim", "fam"):
            p = f"{self.prefix}.{ext}.tmp"
            if os.path.exists(p):
                os.remove(p)


# ---------------------------------------------------------------------------
# Concatenate / subset
# ---------------------------------------------------------------------------

def _read_body(path: str) -> bytes:
    with open(path, "rb") as fh:
        magic = fh.read(3)
        if magic != MAGIC:
            raise ValueError(f"{path} is not a PLINK .bed (bad magic)")
        return fh.read()


def concat_beds(prefixes: Sequence[str], out_prefix: str) -> None:
    """Concatenate variant rows of several beds that share the same sample order.

    The output ``.fam`` is taken from the first input, so all inputs must have an
    identical sample order (enforced by the caller).
    """
    if not prefixes:
        raise ValueError("concat_beds: no input prefixes")
    fam_lines = read_fam(prefixes[0] + ".fam")
    os.makedirs(os.path.dirname(os.path.abspath(out_prefix)) or ".", exist_ok=True)
    with open(out_prefix + ".bed.tmp", "wb") as bed, open(
        out_prefix + ".bim.tmp", "w", encoding="utf-8"
    ) as bim:
        bed.write(MAGIC)
        for prefix in prefixes:
            with open(prefix + ".bim", encoding="utf-8") as fh:
                for line in fh:
                    if line.strip():
                        bim.write(line if line.endswith("\n") else line + "\n")
            bed.write(_read_body(prefix + ".bed"))
    with open(out_prefix + ".fam.tmp", "w", encoding="utf-8") as fh:
        for line in fam_lines:
            fh.write(line + "\n")
    for ext in ("bed", "bim", "fam"):
        os.replace(f"{out_prefix}.{ext}.tmp", f"{out_prefix}.{ext}")


def subset_bed(src_prefix: str, wanted_ids: Iterable[str], out_prefix: str) -> int:
    """Write the variants of ``src`` whose id is in ``wanted_ids``.

    Returns the number of variants written (0 is valid -> an empty bed).
    """
    wanted = set(wanted_ids)
    fam_lines = read_fam(src_prefix + ".fam")
    n = len(fam_lines)
    stride = bytes_per_variant(n)

    bim_lines = read_bim(src_prefix + ".bim")
    keep_idx = [i for i, line in enumerate(bim_lines) if line.split()[1] in wanted]

    os.makedirs(os.path.dirname(os.path.abspath(out_prefix)) or ".", exist_ok=True)
    with open(src_prefix + ".bed", "rb") as bed_in, open(
        out_prefix + ".bed.tmp", "wb"
    ) as bed_out, open(out_prefix + ".bim.tmp", "w", encoding="utf-8") as bim_out:
        bed_out.write(MAGIC)
        for i in keep_idx:
            bim_out.write(bim_lines[i] + "\n")
            bed_in.seek(3 + i * stride)
            bed_out.write(bed_in.read(stride))
    with open(out_prefix + ".fam.tmp", "w", encoding="utf-8") as fh:
        for line in fam_lines:
            fh.write(line + "\n")
    for ext in ("bed", "bim", "fam"):
        os.replace(f"{out_prefix}.{ext}.tmp", f"{out_prefix}.{ext}")
    return len(keep_idx)
