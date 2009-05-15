#!/bin/bash
#
# Fix localversion.

echo "checking localversion"
branch=`git branch 2>&1 | egrep '^\*' | sed -e 's/^\* //' | sed -e 's/[()]//g' | sed -e 's/ /_/g'`
localversion=`git log --pretty=oneline -1 2>&1 | awk '{print $1}'`
if test -n "${branch}" -a -n "${localversion}"; then
  echo "branch=${branch}" > localversion.temp
  echo "localversion=${localversion}" >> localversion.temp
  if test -f localversion; then
    if test -z "`diff localversion localversion.temp`"; then
      echo "localversion did not change"
      rm localversion.temp
    else
      echo "localversion updated"
      diff -u localversion localversion.temp
      mv -f localversion.temp localversion
    fi
  else
    echo "localversion did not exist"
    mv -f localversion.temp localversion
  fi
    rm -f localversion.temp
  echo "localversion at ${branch}:${localversion}"
fi
