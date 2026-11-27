"""Python binding for the mesh-connect shared library (C ABI, via ctypes).

    from meshconnect import Mesh
    m = Mesh("model.bdf")
    welds = m.connect([(x, y, z), ...], tol=1.2)
    print(m.components(), m.components(welds))

The library is found through $MC_LIB, or next to this file / in ../build.
No compiled extension is needed: the same .so / .dll works for any Python.
"""
import ctypes
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))


def _find_lib():
    if os.environ.get("MC_LIB"):
        return os.environ["MC_LIB"]
    names = ["meshconnect.dll"] if sys.platform == "win32" else ["libmeshconnect.so", "libmeshconnect.dylib"]
    dirs = [_HERE, os.path.join(_HERE, "..", "build"), os.path.join(_HERE, "..", "build", "Release")]
    for d in dirs:
        for n in names:
            p = os.path.join(d, n)
            if os.path.exists(p):
                return p
    raise OSError("mesh-connect library not found; set MC_LIB")


class _Weld(ctypes.Structure):
    _fields_ = [
        ("status", ctypes.c_int32),
        ("nsheets", ctypes.c_int32),
        ("elem_id", ctypes.c_int64 * 3),
        ("pid", ctypes.c_int32 * 3),
        ("dist", ctypes.c_double * 3),
        ("point", (ctypes.c_double * 3) * 3),
        ("weight", (ctypes.c_double * 4) * 3),
    ]


_lib = ctypes.CDLL(_find_lib())
_lib.mc_version.restype = ctypes.c_char_p
_lib.mc_last_error.restype = ctypes.c_char_p
_lib.mc_load_nastran.argtypes = [ctypes.c_char_p]
_lib.mc_load_nastran.restype = ctypes.c_void_p
_lib.mc_free.argtypes = [ctypes.c_void_p]
_lib.mc_node_count.argtypes = [ctypes.c_void_p]
_lib.mc_node_count.restype = ctypes.c_int64
_lib.mc_elem_count.argtypes = [ctypes.c_void_p]
_lib.mc_elem_count.restype = ctypes.c_int64
_lib.mc_components.argtypes = [ctypes.c_void_p]
_lib.mc_components.restype = ctypes.c_int32
_lib.mc_connect.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_int64, ctypes.c_double,
                            ctypes.c_int32, ctypes.c_int32, ctypes.POINTER(_Weld)]
_lib.mc_connect.restype = ctypes.c_int64
_lib.mc_components_with.argtypes = [ctypes.c_void_p, ctypes.POINTER(_Weld), ctypes.c_int64]
_lib.mc_components_with.restype = ctypes.c_int32
_lib.mc_write_cwelds.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.POINTER(ctypes.c_double),
                                 ctypes.POINTER(_Weld), ctypes.c_int64, ctypes.c_int64, ctypes.c_int32]
_lib.mc_write_cwelds.restype = ctypes.c_int64

STATUS = {0: "ok", 1: "one_sheet", 2: "no_sheet"}


def version():
    return _lib.mc_version().decode()


class Welds:
    """Result of Mesh.connect: keeps the raw C array for later calls."""

    def __init__(self, xyz, arr):
        self._xyz, self._arr = xyz, arr

    def __len__(self):
        return len(self._arr)

    def __getitem__(self, i):
        w = self._arr[i]
        n = w.nsheets
        return {
            "status": STATUS[w.status],
            "sheets": [{"pid": w.pid[s], "elem": w.elem_id[s], "dist": w.dist[s],
                        "point": tuple(w.point[s]), "weights": tuple(w.weight[s])} for s in range(n)],
        }


class Mesh:
    def __init__(self, path):
        self._h = _lib.mc_load_nastran(os.fsencode(path))
        if not self._h:
            raise RuntimeError(_lib.mc_last_error().decode())

    def close(self):
        if self._h:
            _lib.mc_free(self._h)
            self._h = None

    __del__ = close

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    @property
    def nodes(self):
        return _lib.mc_node_count(self._h)

    @property
    def elems(self):
        return _lib.mc_elem_count(self._h)

    def components(self, welds=None):
        if welds is None:
            return _lib.mc_components(self._h)
        return _lib.mc_components_with(self._h, welds._arr, len(welds))

    def connect(self, points, tol=1.2, max_sheets=3, threads=0):
        n = len(points)
        xyz = (ctypes.c_double * (3 * n))(*[c for p in points for c in p])
        arr = (_Weld * n)()
        if _lib.mc_connect(self._h, xyz, n, tol, max_sheets, threads, arr) < 0:
            raise RuntimeError(_lib.mc_last_error().decode())
        return Welds(xyz, arr)

    def write_cwelds(self, path, welds, first_id, pwid=9000):
        n = _lib.mc_write_cwelds(self._h, os.fsencode(path), welds._xyz, welds._arr, len(welds), first_id, pwid)
        if n < 0:
            raise RuntimeError(_lib.mc_last_error().decode())
        return n
