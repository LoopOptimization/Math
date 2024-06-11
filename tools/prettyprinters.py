from __future__ import print_function

import gdb.printing


class VectorPrinter(Iterator):
    """Print a Vector<>"""

    def __init__(self, val):
        self.val = val
        t = val.type.template_argument(0).pointer()
        self.begin = val["ptr"].cast(t)
        self.size = val["sz"]["value_"].cast(gdb.lookup_type('long'))
        self.i = 0

    def __next__(self):
        if self.i == self.size:
            raise StopIteration
        ret = "[{}]".format(self.i), (self.begin + self.i).dereference()
        self.i += 1
        return ret

    def to_string(self):
        return "Vector of size: {}".format(self.size)

    def display_hint(self):
        return "array"


class BaseMatrixPrinter(Iterator):
    """Print a StridedMatrix<>"""

    def __init__(self, begin, rows, cols, stride):
        self.begin = begin
        self.rows = rows
        self.cols = cols
        self.stride = stride
        self.r = 0
        self.c = -1

    def __next__(self):
        if (self.rows == 0) or (self.cols == 0):
            raise StopIteration
        self.c += 1
        if self.c == self.cols:
            self.c = 0
            self.r += 1
            if self.r == self.rows:
                raise StopIteration
        ind = "[{}, {}]".format(self.r, self.c)
        r = (self.begin + (self.c + self.r * self.stride)).dereference()
        v = str(r)
        if r >= 0:
            v = " " + v
        if (self.c == 0) and (self.r != 0):
            v = "\n    " + v
        if (self.r == self.rows - 1) and (self.c == self.cols - 1):
            v = v + " "
        return ind, v

    def to_string(self):
        return "Matrix, {} x {}, stride {}:\n".format(self.rows, self.cols, self.stride)

    def display_hint(self):
        return "array"


class SquareMatrixPrinter(BaseMatrixPrinter):
    """Print a Matrix<>"""

    def __init__(self, val):
        t = val.type.template_argument(0).pointer()
        M = val["sz"]["m_"]["value_"].cast(gdb.lookup_type('long'))
        BaseMatrixPrinter.__init__(self, val["ptr"].cast(t), M, M, M)


class DenseMatrixPrinter(BaseMatrixPrinter):
    """Print a Matrix<>"""

    def __init__(self, val):
        t = val.type.template_argument(0).pointer()
        M = val["sz"]["m_"]["value_"].cast(gdb.lookup_type('long'))
        N = val["sz"]["n_"]["value_"].cast(gdb.lookup_type('long'))
        BaseMatrixPrinter.__init__(self, val["ptr"].cast(t), M, N, N)


class StridedMatrixPrinter(BaseMatrixPrinter):
    """Print a Matrix<>"""

    def __init__(self, val):
        t = val.type.template_argument(0).pointer()
        BaseMatrixPrinter.__init__(
            self,
            val["ptr"].cast(t),
            val["sz"]["m_"]["value_"].cast(gdb.lookup_type('long')),
            val["sz"]["n_"]["value_"].cast(gdb.lookup_type('long')),
            val["sz"]["stride_m_"]["value_"].cast(gdb.lookup_type('long')),
        )


class WrappedIntegerPrinter:
    """Print an integer wrapped in a struct."""

    def __init__(self, val):
        self.val = str(val["value_"])

    def to_string(self):
        return self.val


pp = gdb.printing.RegexpCollectionPrettyPrinter("PolyMath")
pp.add_printer(
    "poly::math::PtrVector",
    "^poly::math::(Mut)?Array<.*, poly::math::axis::Length<-1l, long>,[^,]*?>$",
    VectorPrinter
)
pp.add_printer(
    "poly::math::Vector",
    "^poly::math::ManagedArray<.*, poly::math::axis::Length<-1l, long>,[^,]*?,[^,]*?>$",
    VectorPrinter,
)
pp.add_printer(
    "poly::math::SquarePtrMatrix",
    "^poly::math::(Mut)?Array<.*, poly::math::SquareDims<-1l>,[^,]*?>$",
    SquareMatrixPrinter,
)
pp.add_printer(
    "poly::math::DensePtrMatrix",
    "^poly::math::(Mut)?Array<.*, poly::math::DenseDims<-1l>,[^,]*?>$",
    DenseMatrixPrinter,
)
pp.add_printer(
    "poly::math::StridedPtrMatrix",
    "^poly::math::(Mut)?Array<.*, poly::math::StridedDims<-1l>,[^,]*?>$",
    StridedMatrixPrinter,
)
pp.add_printer(
    "poly::math::SquareMatrix",
    "^poly::math::ManagedArray<.*, poly::math::SquareDims<-1l>,[^,]*?,[^,]*?>$",
    SquareMatrixPrinter,
)
pp.add_printer(
    "poly::math::DenseMatrix",
    "^poly::math::ManagedArray<.*, poly::math::DenseDims<-1l>,[^,]*?,[^,]*?>$",
    DenseMatrixPrinter,
)
pp.add_printer(
    "poly::math::StridedMatrix",
    "^poly::math::ManagedArray<.*, poly::math::StridedDims<-1l>,[^,]*?,[^,]*?>$",
    StridedMatrixPrinter,
)
pp.add_printer(
    "poly::math::axis::Length", "poly::math::axis::Length<-1l, long>", WrappedIntegerPrinter
)
pp.add_printer(
    "poly::math::axis::Capacity", "poly::math::axis::Capacity<-1l, long>", WrappedIntegerPrinter
)
pp.add_printer(
    "poly::math::axis::Row", "poly::math::axis::Row<-1l>", WrappedIntegerPrinter
)
pp.add_printer(
    "poly::math::axis::Col", "poly::math::axis::Col<-1l>", WrappedIntegerPrinter
)
pp.add_printer(
    "poly::math::axis::RowStride",
    "poly::math::axis::RowStride<-1l>",
    WrappedIntegerPrinter,
)
pp.add_printer(
    "poly::math::axis::Capacity",
    "poly::math::axis::Capacity<-1l>",
    WrappedIntegerPrinter,
)
gdb.printing.register_pretty_printer(gdb.current_objfile(), pp)
