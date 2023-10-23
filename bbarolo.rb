class Bbarolo < Formula
  desc "3D fitting tool to derive the kinematics of galaxies"
  homepage "https://editeodoro.github.io/Bbarolo/"

  # Default is v1.7 stable
  stable do
    url "https://github.com/editeodoro/Bbarolo/archive/1.7.tar.gz"
    sha256 "445b7b6eae023efca60e669b3e40a1298846772bc2d23a24d69380e52592fa61"
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
  depends_on "gnuplot" => :optional

  def install
    # BBarolo requires a c++11 compiler
    ENV.cxx11

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