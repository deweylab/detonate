mkdir boost
~/downloads/boost/boost_1_52_0/bin.v2/tools/bcp/gcc-4.4.7/release/link-static/bcp --boost=$HOME/downloads/boost/boost_1_52_0 --scan $(ls *.cpp *.hh) boost
~/downloads/boost/boost_1_52_0/bin.v2/tools/bcp/gcc-4.4.7/release/link-static/bcp --boost=$HOME/downloads/boost/boost_1_52_0 build boost
