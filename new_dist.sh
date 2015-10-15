#!/bin/bash

# Default package name...
pack_name=Inelastica

tmp_dir=.tmp
# Get the revision number...
r=`svn info | grep Revision | awk '{print $2}'`
echo "Found revision $r"

function exp_clean {
    local dir=$1 ; shift

    # Clean up the temporary directory
    rm -rf $tmp_dir
    mkdir $tmp_dir
    
    # Export current revision to the directory folder
    svn export . $tmp_dir/$dir
    
    # Remove certain files from the distribution
    for f in new_dist.sh ; do
	rm -r $tmp_dir/$dir/$f
    done
}

function dist_tar {
    local out=$1 ; shift
    local dir=$1 ; shift
    pushd $tmp_dir > /dev/null
    tar cfz ../$out.tar.gz $dir
    popd > /dev/null
}

# Create a new distribution.
# It will ONLY take the latest committed revision.
exp_clean Inelastica-$r
dist_tar $pack_name-$r Inelastica-$r

while [ $# -gt 0 ]; do
    
    # Create all user distributions
    exp_clean Inelastica-$1
    dist_tar $pack_name-$1 Inelastica-$1
    
    shift

done
    
rm -rf $tmp_dir
