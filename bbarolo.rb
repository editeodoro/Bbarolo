class Bbarolo < Formula
  desc "3D fitting tool to derive the kinematics of galaxies"
  homepage "https://editeodoro.github.io/Bbarolo/"

  # Default is v1.7 stable
  stable do
    url "https://github.com/editeodoro/Bbarolo/archive/1.7.tar.gz"
    #sha256 "5d169ffb0e60042c74bf54c6f44ed1d7f3d9a8b871f6c6384a24078c6fe2b172"
  end

  # To install instead latest non-stable version, use --devel option
  head do
    url "https://github.com/editeodoro/Bbarolo/archive/master.tar.gz"
    #sha256 "fb91ff62a33dbda6fb2c9da04e0746017d4e905e981929b8ca7d1336fff0810a"
    version "1.7dev"
  end

  # Dependencies 
  depends_on "cfitsio"
  depends_on "fftw"
  depends_on "wcslib"
  depends_on "gcc"
  depends_on "gnuplot" => :optional
  
  # With Clang 14 gives a segfault, using gcc instead
  fails_with :clang do
    build 1400
    cause "Miscompilation resulting in segfault on queries"
  end
  
  def install
    # BBarolo requires a c++17 compiler
    ENV.cxx17

    # Configure script arguments
    args = %W[
      --prefix=#{prefix}
      --with-cfitsio=#{Formula["cfitsio"].opt_prefix}
      --with-fftw3=#{Formula["fftw"].opt_prefix}
      --with-wcslib=#{Formula["wcslib"].opt_prefix}
    ]

    system "./configure", *args
    system "make"
    system "make", "install"
  end

  test do
    assert_match "This is BBarolo", `#{bin}/BBarolo -v`
  end
end