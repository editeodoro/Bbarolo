class Bbarolo < Formula
  desc "3D fitting tool to derive the kinematics of galaxies"
  homepage "https://editeodoro.github.io/Bbarolo/"

  # Default is v1.6 stable
  stable do
    url "https://github.com/editeodoro/Bbarolo/archive/1.6.tar.gz"
    sha256 "f9ad7d7a32b9141e76d21867ebe17fde777e40a998e08290785d8a284bc3973d"
  end

  # To install instead latest non-stable version, use --devel option
  head do
    url "https://github.com/editeodoro/Bbarolo/archive/master.tar.gz"
    sha256 "3db5d9080b17e2f5b1d70954db2fa06191c9fd9cac31132f0dc48019d7a9f7e2"
    version "1.6dev"
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