class Bbarolo < Formula
  desc "3D fitting tool to derive the kinematics of galaxies"
  homepage "https://editeodoro.github.io/Bbarolo/"
  
  # Default is v1.5 stable
  stable do
    url "https://github.com/editeodoro/Bbarolo/archive/1.5.tar.gz"
    sha256 "38276adf408de83b8aa72e52f022709c3955e98fa378e9d414a2b82e59c4a3c2"
  end
  
  # To install instead latest non-stable version, use --devel option
  devel do
    url "https://github.com/editeodoro/Bbarolo/archive/master.tar.gz"
    sha256 "2d2c8c9a816cf9ba30d2005e8c54141ac5f42461be7b91fdc9f89fd4cc638488"
    version "1.5dev"
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
