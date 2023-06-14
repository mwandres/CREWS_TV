ADCIRC+SWAN installation CygWin

1) download and install cygwin with wget only (download cygwin https://www.cygwin.com/, and during installation select install wget option)

2) wget https://raw.githubusercontent.com/transcode-open/apt-cyg/master/apt-cyg 

chmod +x apt-cyg
mv apt-cyg /usr/bin/

3) Install Dependencies

apt-cyg install _autorebase alternatives appstream at-spi2-core base-cygwin base-files bash binutils bzip2 ca-certificates cmake cmake-gui coreutils crypto-policies cygutils cygwin cygwin-devel dash dconf-service dejavu-fonts desktop-file-utils diffutils dri-drivers editrights extra-cmake-modules file findutils gamin gawk gcc-core gcc-fortran gcc-g++ gdk-pixbuf2.0-svg getent glib2.0-networking gnupg grep groff gsettings-desktop-schemas gtk-update-icon-cache gzip hicolor-icon-theme hostname info ipc-utils less libGL1 libICE6 libQt5Core5 libQt5Gui5 libSM6 libX11-xcb1 libX11_6 libXau6 libXcomposite1 libXcursor1 libXdamage1 libXdmcp6 libXext6 libXfixes3 libXft2 libXi6 libXinerama1 libXrandr2 libXrender1 libXtst6 libappstream4 libarchive13 libargp libassuan0 libatk-bridge2.0_0 libatk1.0_0 libatomic1 libatspi0 libattr1 libblkid1 libbrotlicommon1 libbrotlidec1 libbz2_1 libcairo2 libcares2 libclang8 libcom_err2 libcroco0.6_3 libcrypt2 libcurl4 libdatrie1 libdb5.3 libdbus1_3 libdeflate0 libedit0 libepoxy0 libevent2.0_5 libexempi-devel libexempi3 libexpat1 libfam0 libfdisk1 libffi6 libfontconfig-common libfontconfig1 libfreetype6 libgc1 libgcc1 libgdbm6 libgdbm_compat4 libgdk_pixbuf2.0_0 libgfortran4 libgfortran5 libglapi0 libglib2.0_0 libgmp10 libgnutls30 libgomp1 libgpg-error0 libgpgme11 libgraphite2_3 libgssapi_krb5_2 libgtk3_0 libguile2.2_1 libharfbuzz0 libhdf5_101 libhdf5_103 libhdf5hl_100 libhogweed4 libhwloc15 libiconv2 libicu61 libidn2_0 libintl8 libisl15 libisl22 libjasper4 libjbig2 libjpeg8 libjson-glib1.0_0 libjsoncpp24 libk5crypto3 libkrb5_3 libkrb5support0 libllvm8 libltdl7 liblz4_1 liblzma5 liblzo2_2 libmetalink3 libmpc3 libmpfr6 libncursesw10 libnetcdf-devel libnetcdf-fortran-devel libnetcdf-fortran_6 libnetcdf-fortran_7 libnetcdf13 libnetcdf15 libnettle6 libnghttp2_14 libopenldap2_4_2 libopenmpi-devel libopenmpi12 libopenmpi40 libopenmpicxx1 libopenmpifh12 libopenmpifh40 libopenmpiusef08_40 libopenmpiusetkr40 libp11-kit0 libpango1.0_0 libpcre1 libpcre2_16_0 libpcre2_8_0 libpipeline1 libpixman1_0 libpkgconf3 libpng16 libpolly8 libpopt-common libpopt0 libproxy1 libpsl5 libquadmath0 libreadline7 librest0.7_0 librhash0 librsvg2_2 libsasl2_3 libsigsegv2 libsmartcols1 libsoup-gnome2.4_1 libsoup2.4_1 libsqlite3_0 libssh2_1 libssl1.0 libssl1.1 libstdc++6 libtasn1_6 libthai0 libtiff6 libunistring2 libusb0 libuuid-devel libuuid1 libuv1 libwebp7 libxcb-glx0 libxcb-icccm4 libxcb-image0 libxcb-keysyms1 libxcb-randr0 libxcb-render-util0 libxcb-render0 libxcb-shape0 libxcb-shm0 libxcb-sync1 libxcb-util1 libxcb-xfixes0 libxcb-xinerama0 libxcb-xkb1 libxcb1 libxkbcommon0 libxml2 libyaml0_2 libzstd1 login make man-db mintty ncurses netcdf openmpi openssl p11-kit p11-kit-trust perl perl_autorebase perl_base pkg-config pkgconf publicsuffix-list-dafsa python2 python27 python27-clang python27-pip python27-setuptools python38 python38-pip python38-setuptools rebase run sed shared-mime-info tar terminfo terminfo-extra tzcode tzdata util-linux vim vim-common vim-minimal w32api-headers w32api-runtime wget which windows-default-manifest xkeyboard-config xxd xz zlib0 zstd


4) Extract Source

tar xvf  adcirc_v53_dev_050217.tar
mkdir build
cd build

5) Compile

   option 1: cmake .. -DBUILD_ADCIRC=ON -DBUILD_PADCIRC=ON -DBUILD_ADCSWAN=ON -DBUILD_PADCSWAN=ON -DBUILD_ADCPREP=ON -DBUILD_UTILITIES=OFF -DBUILD_ASWIP=ON -DBUILD_SWAN=ON -DBUILD_PUNSWAN=ON -DENABLE_OUTPUT_NETCDF=ON -DENABLE_OUTPUT_XDMF=OFF -DNETCDFHOME=/usr -DBUILD_UTILITIES=ON -DCMAKE_Fortran_FLAGS="-mtune=native -std=legacy"
             make

   option 2: ccmake 


Anaconda 3 installation

1) Download and install anaconda3: https://docs.anaconda.com/anaconda/install/hashes/win-3-64/
2) Create virtual environment: conda create -n Tuvalu_EWS python=3.9.4 spyder=5.0.3
3) conda activate Tuvalu_EWS

4) install modules

	conda install -c conda-forge pytmd==1.0.4
        conda update pytmd
	conda install -c conda-forge fiona
        pip install adcircpy==1.0.20
	conda install -c conda-forge motuclient
	conda install -c conda-forge geopy
	conda install scikit-learn
	conda istall xarray

5) Install Java

Run Tarawa_EWS

1) Double click Run_operational.bat
2) Windows task scheduler calls Run_operational.bat every 6 hours





