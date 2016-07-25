#!/usr/bin/env python

#usage: make -j 4 2> MAKE ; cat MAKE | python parse_gcc_output.py | less -R
import sys,fileinput

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

import re

for line in fileinput.input():
    line = line.replace("gatb::core::debruijn::impl::Edge_t<gatb::core::debruijn::impl::Node_t<gatb::core::tools::math::IntegerTemplate<boost::mpl::vector4<mpl_::int_<32>, mpl_::int_<64>, mpl_::int_<96>, mpl_::int_<128> > > > >", bcolors.OKBLUE + "Edge" + bcolors.ENDC) # default KSIZE_LIST
    line = line.replace("gatb::core::debruijn::impl::Edge_t<gatb::core::debruijn::impl::Node_t<gatb::core::tools::math::IntegerTemplate<boost::mpl::vector1<mpl_::int_<32> > > > >", bcolors.OKBLUE + "Edge" + bcolors.ENDC) # when compiled with KSIZE_LIST=32
    line = line.replace("gatb::core::debruijn::impl::Edge_t<gatb::core::debruijn::impl::Node_t<gatb::core::tools::math::LargeInt<1> > >", bcolors.OKBLUE + "EdgeFast<1>" + bcolors.ENDC)

    line = line.replace("gatb::core::debruijn::impl::Node_t<gatb::core::tools::math::IntegerTemplate<boost::mpl::vector4<mpl_::int_<32>, mpl_::int_<64>, mpl_::int_<96>, mpl_::int_<128> > > >", bcolors.OKBLUE + "Node" + bcolors.ENDC)
    line = line.replace("gatb::core::debruijn::impl::Node_t<gatb::core::tools::math::IntegerTemplate<boost::mpl::vector1<mpl_::int_<32> > > >", bcolors.OKBLUE + "Node" + bcolors.ENDC)
    line = line.replace("gatb::core::debruijn::impl::Node_t<gatb::core::tools::math::LargeInt<1> >", bcolors.OKBLUE + "NodeFast<1>" + bcolors.ENDC)

    line = line.replace("boost::variant<boost::detail::variant::over_sequence<boost::mpl::v_item<gatb::core::debruijn::impl::GraphData<128ul>, boost::mpl::v_item<gatb::core::debruijn::impl::GraphData<96ul>, boost::mpl::v_item<gatb::core::debruijn::impl::GraphData<64ul>, boost::mpl::v_item<gatb::core::debruijn::impl::GraphData<32ul>, boost::mpl::vector0<mpl_::na>, 0>, 0>, 0>, 0> >> ", bcolors.OKBLUE + "GraphDataVariant" + bcolors.ENDC)
    line = line.replace("boost::variant<boost::detail::variant::over_sequence<boost::mpl::v_item<gatb::core::debruijn::impl::GraphData<32ul>, boost::mpl::vector0<mpl_::na>, 0> >>", bcolors.OKBLUE + "GraphDataVariant" + bcolors.ENDC)
    line = line.replace("boost::variant<gatb::core::debruijn::impl::GraphData<32ul>>",  bcolors.OKBLUE + "GraphDataVariantFast<1>" + bcolors.ENDC)
    line = line.replace("boost::variant<gatb::core::debruijn::impl::GraphData<32ul> >", bcolors.OKBLUE + "GraphDataVariantFast<1>" + bcolors.ENDC)

    line = line.replace("boost::variant<boost::detail::variant::over_sequence<boost::mpl::v_item<gatb::core::debruijn::impl::GraphData<128ul>, boost::mpl::v_item<gatb::core::debruijn::impl::GraphData<96ul>, boost::mpl::v_item<gatb::core::debruijn::impl::GraphData<64ul>, boost::mpl::v_item<gatb::core::debruijn::impl::GraphData<32ul>, boost::mpl::vector0<mpl_::na>, 0>, 0>, 0>, 0> > >", bcolors.OKBLUE + "GraphDataVariant" + bcolors.ENDC)
    line = line.replace("boost::variant<boost::detail::variant::over_sequence<boost::mpl::v_item<gatb::core::debruijn::impl::GraphData<128>, boost::mpl::v_item<gatb::core::debruijn::impl::GraphData<96>, boost::mpl::v_item<gatb::core::debruijn::impl::GraphData<64>, boost::mpl::v_item<gatb::core::debruijn::impl::GraphData<32>, boost::mpl::vector0<mpl_::na>, 0>, 0>, 0>, 0> >",  bcolors.OKBLUE + "GraphDataVariant" + bcolors.ENDC) # clang flavor

    line = line.replace("gatb::core::tools::math::LargeInt<1>", bcolors.OKBLUE + "LargeInt<1>" + bcolors.ENDC)
    line = line.replace("gatb::core::debruijn::impl::GraphTemplate", bcolors.OKBLUE + "GraphTemplate" + bcolors.ENDC)
    line = line.replace("gatb::core::debruijn::impl::GraphUnitigsTemplate", bcolors.OKBLUE + "GraphUnitigsTemplate" + bcolors.ENDC)

    line = line.replace("undefined reference to ", bcolors.FAIL + "undefined reference to " + bcolors.ENDC)
    line = line.replace("In function", bcolors.WARNING + "In function" + bcolors.ENDC)
    line = line.replace("gatb::core::tools::math::IntegerTemplate<boost::mpl::vector4<mpl_::int_<32>, mpl_::int_<64>, mpl_::int_<96>, mpl_::int_<128> > > ", bcolors.OKBLUE + "Integer" + bcolors.ENDC)
    line = re.sub(r"([^\:]*):\(.text.[^\)]*\)",  bcolors.OKGREEN + r"\1:(..)" + bcolors.ENDC,line)

    sys.stdout.write(line)

