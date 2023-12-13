#run rSFS
NITER=10
echo "running mstatpop option 92 (rSFS values) for DIPLOIDS with 4 different conditions and $NITER iterations ..."
./mstatspop -f ms -i ./slim.snm.output_all.ms -o 92 -N 50 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 -n ./chr_ext.txt -l 500000 -r $NITER -F 1 -T ./slim.snm.output_file.ms_mstatspopo92.txt
./mstatspop -f ms -i ./slim.selsweep.output_all.ms -o 92 -N 50 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 -n ./chr_ext.txt -l 500000 -r $NITER -F 1 -T ./slim.selsweep.output_file.ms_mstatspopo92.txt
./mstatspop -f ms -i ./slim.incomplete_selsweep.output_all.ms -o 92 -N 50 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 -n ./chr_ext.txt -l 500000 -r $NITER -F 1 -T ./slim.incomplete_selsweep.output_file.ms_mstatspopo92.txt
./mstatspop -f ms -i ./slim.standing_selsweep.output_all.ms -o 92 -N 50 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 -n ./chr_ext.txt -l 500000 -r $NITER -F 1 -T ./slim.standing_selsweep.output_file.ms_mstatspopo92.txt

echo "running mstatpop option 92 (rSFS values) for HAPLOIDS with 4 different conditions and $NITER iterations ..."
./mstatspop -f ms -i ./slim.snm.output_all.ms -o 92 -N 100 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 -n ./chr_ext.txt -l 500000 -r $NITER -F 1 -T ./slim.snm.output_file.ms_hap_mstatspopo92.txt
./mstatspop -f ms -i ./slim.selsweep.output_all.ms -o 92 -N 100 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 -n ./chr_ext.txt -l 500000 -r $NITER -F 1 -T ./slim.selsweep.output_file.ms_hap_mstatspopo92.txt
./mstatspop -f ms -i ./slim.incomplete_selsweep.output_all.ms -o 92 -N 100 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 -n ./chr_ext.txt -l 500000 -r $NITER -F 1 -T ./slim.incomplete_selsweep.output_file.ms_hap_mstatspopo92.txt
./mstatspop -f ms -i ./slim.standing_selsweep.output_all.ms -o 92 -N 100 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 -n ./chr_ext.txt -l 500000 -r $NITER -F 1 -T ./slim.standing_selsweep.output_file.ms_hap_mstatspopo92.txt

echo "running mstatpop option 1 for WHOLE POPULATION with 4 different conditions and $NITER iterations ..."
#run stats
./mstatspop -f ms -i ./slim.snm.output_all.ms -o 1 -N 1 100 -n ./chr_ext.txt -l 500000 -r $NITER -F 1 -T ./slim.snm.output_file.ms_mstatspopo1.txt
./mstatspop -f ms -i ./slim.selsweep.output_all.ms -o 1 -N 1 100 -n ./chr_ext.txt -l 500000 -r $NITER -F 1 -T ./slim.selsweep.output_file.ms_mstatspopo1.txt
./mstatspop -f ms -i ./slim.incomplete_selsweep.output_all.ms -o 1 -N 1 100 -n ./chr_ext.txt -l 500000 -r $NITER -F 1 -T ./slim.incomplete_selsweep.output_file.ms_mstatspopo1.txt
./mstatspop -f ms -i ./slim.standing_selsweep.output_all.ms -o 1 -N 1 100 -n ./chr_ext.txt -l 500000 -r $NITER -F 1 -T ./slim.standing_selsweep.output_file.ms_mstatspopo1.txt

echo "keep desired columns with 4 different conditions and $NITER iterations ..."
#keep desired columns
perl ./collect_data_columns.pl -in  ./slim.snm.output_file.ms_mstatspopo1.txt  -fc columns_to_choose.txt > ./slim.snm.output_file.ms_mstatspopo1.txt.columns.txt
perl ./collect_data_columns.pl -in  ./slim.selsweep.output_file.ms_mstatspopo1.txt   -fc columns_to_choose.txt > ./slim.selsweep.output_file.ms_mstatspopo1.txt.columns.txt
perl ./collect_data_columns.pl -in  ./slim.incomplete_selsweep.output_file.ms_mstatspopo1.txt   -fc columns_to_choose.txt > ./slim.incomplete_selsweep.output_file.ms_mstatspopo1.txt.columns.txt
perl ./collect_data_columns.pl -in  ./slim.standing_selsweep.output_file.ms_mstatspopo1.txt  -fc columns_to_choose.txt > ./slim.standing_selsweep.output_file.ms_mstatspopo1.txt.columns.txt

-----------------------------------------------------------------------
#run rSFS
NITER=100
echo "running mstatpop option 92 (rSFS values) for DIPLOIDS with 4 different conditions and $NITER iterations ..."
./mstatspop -f ms -i ./slim.snm.output_all.ms -o 92 -N 50 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 -n ./chr_ext.txt -l 100000 -r $NITER -F 1 -T ./slim.snm.output_file.ms_mstatspopo92.txt
./mstatspop -f ms -i ./slim.selsweep.output_all.ms -o 92 -N 50 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 -n ./chr_ext.txt -l 100000 -r $NITER -F 1 -T ./slim.selsweep.output_file.ms_mstatspopo92.txt
./mstatspop -f ms -i ./slim.incomplete_selsweep.output_all.ms -o 92 -N 50 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 -n ./chr_ext.txt -l 100000 -r $NITER -F 1 -T ./slim.incomplete_selsweep.output_file.ms_mstatspopo92.txt
./mstatspop -f ms -i ./slim.standing_selsweep.output_all.ms -o 92 -N 50 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 -n ./chr_ext.txt -l 100000 -r $NITER -F 1 -T ./slim.standing_selsweep.output_file.ms_mstatspopo92.txt

echo "running mstatpop option 92 (rSFS values) for HAPLOIDS with 4 different conditions and $NITER iterations ..."
./mstatspop -f ms -i ./slim.snm.output_all.ms -o 92 -N 100 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 -n ./chr_ext.txt -l 500000 -r $NITER -F 1 -T ./slim.snm.output_file.ms_hap_mstatspopo92.txt
./mstatspop -f ms -i ./slim.selsweep.output_all.ms -o 92 -N 100 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 -n ./chr_ext.txt -l 500000 -r $NITER -F 1 -T ./slim.selsweep.output_file.ms_hap_mstatspopo92.txt
./mstatspop -f ms -i ./slim.incomplete_selsweep.output_all.ms -o 92 -N 100 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 -n ./chr_ext.txt -l 500000 -r $NITER -F 1 -T ./slim.incomplete_selsweep.output_file.ms_hap_mstatspopo92.txt
./mstatspop -f ms -i ./slim.standing_selsweep.output_all.ms -o 92 -N 100 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 -n ./chr_ext.txt -l 500000 -r $NITER -F 1 -T ./slim.standing_selsweep.output_file.ms_hap_mstatspopo92.txt

echo "running mstatpop option 1 for WHOLE POPULATION with 4 different conditions and $NITER iterations ..."
#run stats
./mstatspop -f ms -i ./slim.snm.output_all.ms -o 1 -N 1 100 -n ./chr_ext.txt -l 100000 -r $NITER -F 1 -T ./slim.snm.output_file.ms_mstatspopo1.txt
./mstatspop -f ms -i ./slim.selsweep.output_all.ms -o 1 -N 1 100 -n ./chr_ext.txt -l 100000 -r $NITER -F 1 -T ./slim.selsweep.output_file.ms_mstatspopo1.txt
./mstatspop -f ms -i ./slim.incomplete_selsweep.output_all.ms -o 1 -N 1 100 -n ./chr_ext.txt -l 100000 -r $NITER -F 1 -T ./slim.incomplete_selsweep.output_file.ms_mstatspopo1.txt
./mstatspop -f ms -i ./slim.standing_selsweep.output_all.ms -o 1 -N 1 100 -n ./chr_ext.txt -l 100000 -r $NITER -F 1 -T ./slim.standing_selsweep.output_file.ms_mstatspopo1.txt

echo "keep desired columns with 4 different conditions and $NITER iterations ..."
#keep desired columns
perl ./collect_data_columns.pl -in  ./slim.snm.output_file.ms_mstatspopo1.txt  -fc columns_to_choose.txt > ./slim.snm.output_file.ms_mstatspopo1.txt.columns.txt
perl ./collect_data_columns.pl -in  ./slim.selsweep.output_file.ms_mstatspopo1.txt   -fc columns_to_choose.txt > ./slim.selsweep.output_file.ms_mstatspopo1.txt.columns.txt
perl ./collect_data_columns.pl -in  ./slim.incomplete_selsweep.output_file.ms_mstatspopo1.txt   -fc columns_to_choose.txt > ./slim.incomplete_selsweep.output_file.ms_mstatspopo1.txt.columns.txt
perl ./collect_data_columns.pl -in  ./slim.standing_selsweep.output_file.ms_mstatspopo1.txt  -fc columns_to_choose.txt > ./slim.standing_selsweep.output_file.ms_mstatspopo1.txt.columns.txt
