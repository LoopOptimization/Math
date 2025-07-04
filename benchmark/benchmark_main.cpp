import Nanobench;
import std;

// Dual benchmark functions
template <std::ptrdiff_t M, std::ptrdiff_t N> void BM_dualprod(Bench &bench);

template <std::ptrdiff_t M, std::ptrdiff_t N>
void BM_dualprod_manual(Bench &bench);

template <std::ptrdiff_t M, std::ptrdiff_t N>
void BM_dualprod_simdarray(Bench &bench);

template <std::ptrdiff_t M, std::ptrdiff_t N>
void BM_dualprod_manual_tuple(Bench &bench);

template <std::ptrdiff_t M, std::ptrdiff_t N>
void BM_dualprod_simdarray_tuple(Bench &bench);

template <std::ptrdiff_t M, std::ptrdiff_t N>
void BM_dualdivsum(Bench &bench, std::ptrdiff_t len);

// Simplex benchmark functions
void BM_Simplex0(Bench &bench);
void BM_Simplex1(Bench &bench);

int main() {
  Bench bench;

  // Dual product benchmarks
  bench.title("Dual Product Benchmarks");
  bench.run("BM_dualprod<6,2>", [&] { BM_dualprod<6, 2>(bench); });
  bench.run("BM_dualprod<7,2>", [&] { BM_dualprod<7, 2>(bench); });
  bench.run("BM_dualprod<8,2>", [&] { BM_dualprod<8, 2>(bench); });

  bench.run("BM_dualprod_manual<6,2>",
            [&] { BM_dualprod_manual<6, 2>(bench); });
  bench.run("BM_dualprod_manual<7,2>",
            [&] { BM_dualprod_manual<7, 2>(bench); });
  bench.run("BM_dualprod_manual<8,2>",
            [&] { BM_dualprod_manual<8, 2>(bench); });

  bench.run("BM_dualprod_simdarray<6,2>",
            [&] { BM_dualprod_simdarray<6, 2>(bench); });
  bench.run("BM_dualprod_simdarray<7,2>",
            [&] { BM_dualprod_simdarray<7, 2>(bench); });
  bench.run("BM_dualprod_simdarray<8,2>",
            [&] { BM_dualprod_simdarray<8, 2>(bench); });

  bench.run("BM_dualprod_manual_tuple<6,2>",
            [&] { BM_dualprod_manual_tuple<6, 2>(bench); });
  bench.run("BM_dualprod_manual_tuple<7,2>",
            [&] { BM_dualprod_manual_tuple<7, 2>(bench); });
  bench.run("BM_dualprod_manual_tuple<8,2>",
            [&] { BM_dualprod_manual_tuple<8, 2>(bench); });

  bench.run("BM_dualprod_simdarray_tuple<6,2>",
            [&] { BM_dualprod_simdarray_tuple<6, 2>(bench); });
  bench.run("BM_dualprod_simdarray_tuple<7,2>",
            [&] { BM_dualprod_simdarray_tuple<7, 2>(bench); });
  bench.run("BM_dualprod_simdarray_tuple<8,2>",
            [&] { BM_dualprod_simdarray_tuple<8, 2>(bench); });

  // Dual division sum benchmarks with different sizes
  bench.title("Dual Division Sum Benchmarks");
  for (std::ptrdiff_t len = 1; len <= 1024; len *= 2) {
    bench.run("BM_dualdivsum<7,0>_len=" + std::to_string(len),
              [&] { BM_dualdivsum<7, 0>(bench, len); });
    bench.run("BM_dualdivsum<8,0>_len=" + std::to_string(len),
              [&] { BM_dualdivsum<8, 0>(bench, len); });
    bench.run("BM_dualdivsum<7,2>_len=" + std::to_string(len),
              [&] { BM_dualdivsum<7, 2>(bench, len); });
    bench.run("BM_dualdivsum<8,2>_len=" + std::to_string(len),
              [&] { BM_dualdivsum<8, 2>(bench, len); });
    bench.run("BM_dualdivsum<7,4>_len=" + std::to_string(len),
              [&] { BM_dualdivsum<7, 4>(bench, len); });
    bench.run("BM_dualdivsum<8,4>_len=" + std::to_string(len),
              [&] { BM_dualdivsum<8, 4>(bench, len); });
  }

  // Simplex benchmarks
  bench.title("Simplex Benchmarks");
  bench.run("BM_Simplex0", [&] { BM_Simplex0(bench); });
  bench.run("BM_Simplex1", [&] { BM_Simplex1(bench); });

  return 0;
}
