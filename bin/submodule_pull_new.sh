for dir in */; do
  branch="${dir%/}"  # remove trailing slash to get the branch name
  if [ -f "$dir/.git" ]; then
    echo "Processing $dir -> branch $branch"
    (
      cd "$dir" || exit
      git status
      git pull
    )
  fi
done
