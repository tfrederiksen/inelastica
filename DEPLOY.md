
Deployment of Inelastica
========================

This document describes the deployment details to perform a version release.

Version release
---------------

The release cycle should be performed like this:

1. Increment the release numbers in the top-directory `setup.py` script:
   
        MAJOR
        MINOR
        MICRO

   The complete tag is `VERSION = "v[MAJOR].[MINOR].[MICRO]"`.
	
2. Set the variable `GIT_REVISION` to the latest commit in the top-directory `setup.py`.
   This means that the revision specification for the release
   actually corresponds to the commit just before the actual release.
   You can get the commit hash by:

        git rev-parse HEAD
        GIT_REVISION = <git rev-parse HEAD>

3. Add `setup.py` to the commit:

        git commit [-s] -m "Inelastica release: [VERSION]"

   with the corresponding version number.

4. Tag the commit with:

        git tag -a "[VERSION]" -m "Releasing [VERSION]"

   Include a complete description of the essential code modifications since last release.

5. Merge `master` branch into `gh-pages` and update `html` documentation:

        git checkout gh-pages
        git merge master
        cd docs
        ./run.sh
        git add latest
        git commit -m "docs: update html"
