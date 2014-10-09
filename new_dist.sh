#!/bin/bash

# Default package name...
pack_name=Inelastica

tmp_dir=.tmp
# Get the revision number...
r=`svn info | grep Revision | awk '{print $2}'`
echo "Found revision $r"

function exp_clean {
    local dir=$1 ; shift
    svn export . $tmp_dir/$dir
    rm -rf $tmp_dir/$dir/dist
}

function dist_tar {
    local out=$1 ; shift
    local dir=$1 ; shift
    pushd $tmp_dir > /dev/null
    tar cfz ../dist/$out.tar.gz $dir
    popd > /dev/null
}

# Create a new distribution.
# It will ONLY take the latest committed revision.
rm -rf $tmp_dir
mkdir $tmp_dir

#exp_clean Inelastica
#dist_tar $pack_name Inelastica

#exp_clean Inelastica-DEV
#dist_tar $pack_name-DEV Inelastica-DEV

exp_clean Inelastica-$r
dist_tar $pack_name-$r Inelastica-$r

rm -rf $tmp_dir