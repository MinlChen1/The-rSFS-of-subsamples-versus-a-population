echo "run slim simulations for SNM and for selective sweeps: Complete selective sweep, incomplete selective sweep, soft selective sweep:"
echo "conditions of the SLiM simulations are included in the script sweeps_ind.slim: mu=2e-7, r=1e-7, L=5e5, 2Ne=1000"
rm *.ms
rm *.ms.log
NITER=10
for i in $(seq 1 ${NITER})
do
    echo iteration $i/$NITER
    slim -t -m -d "freq_sel_init=0.001" -d "freq_sel_end=0.999" -d "s_beneficial=0.00" -d "ind_sample_size=50" -d "file_output1='./slim.snm.output_all.ms'" ./sweeps_ind.slim 1>>  ./slim.snm.output_all.ms.log
    slim -t -m -d "freq_sel_init=0.001" -d "freq_sel_end=1.000" -d "s_beneficial=0.05" -d "ind_sample_size=50" -d "file_output1='./slim.selsweep.output_all.ms'" ./sweeps_ind.slim 1>>  ./slim.selsweep.output_all.ms.log
    slim -t -m -d "freq_sel_init=0.001" -d "freq_sel_end=0.750" -d "s_beneficial=0.05" -d "ind_sample_size=50" -d "file_output1='./slim.incomplete_selsweep.output_all.ms'" ./sweeps_ind.slim 1>> ./slim.incomplete_selsweep.output_all.ms.log
    slim -t -m -d "freq_sel_init=0.500" -d "freq_sel_end=1.000" -d "s_beneficial=0.05" -d "ind_sample_size=50" -d "file_output1='./slim.standing_selsweep.output_all.ms'" ./sweeps_ind.slim 1>>  ./slim.standing_selsweep.output_all.ms.log
done
#rm tmp_slim_*.txt
