## The following lines give status codes and some logging
echo "Job status:$?" >> outputs/log.txt
echo " $(hostname) Finished [% script %] at $(date), it took $(( SECONDS / 60 )) minutes." >> outputs/log.txt
touch [% finished_file %]
