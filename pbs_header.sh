#
# Version 1.0
# Written by Joshua Colvin
# jcolvin@tgen.org
#
# Function used to log every command that is executed.
# Takes two arguments:
#  - Directory command should be executed in.
#  - Command to execute.
#    Note that it is important to surround entire command with quotes
#    so that pipes and redirects are handled properly.
#
# Example:
# echo_and_exec $PBS_O_WORKDIR "wc -l *.fastq > fastq_line_count.txt"
#
echo_and_exec()
{
  cd $1
  echo ""
  echo "In $1"
  shift
  echo "Running: $*"
  if [ ! "$test" ]
  then
    eval $*
  fi
  status=$?

  if [ ! $status == 0 ]
  then
    echo "Previous command returned error: $status"
    exit 1;
  fi
}

# List of variable names that are required to continue
# (names only, do not include '$')
#
# Example:
# require_var foo bar
#
require_var()
{
  for var in $*
  do
    eval if [ ! "\"\$$var\"" ]\; then echo "Oops, need variable $var\!"\; failed=1\; fi\;
  done

  if [ "$failed" ]
  then
    exit 1;
  fi

}

if [ ! $PBS_O_WORKDIR ]
then
  # Called from command line, not in PBS.
  PBS_O_WORKDIR="$PWD"
fi

cd ${PBS_O_WORKDIR}