#!/bin/sh

set -ex
aclocal -I m4
autoconf
chmod a+x configure
autoheader
automake -a
