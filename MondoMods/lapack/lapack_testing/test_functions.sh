function run_test
{
  $1 < $2 > $3 2>&1

  if test `egrep -i "(fail)" $3 | wc -l` -ne 0; then
    echo "test $1 with input $2 and result in $3 failed"
    return 1
  else
    return 0
  fi
}
