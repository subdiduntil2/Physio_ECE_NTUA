function tumour_GC=GC_out_of_cancer_final(i,j,k,m,time)
 
 max_side_of_mesh=(m-1)/2;
 
tumour_GC.coords=[-max_side_of_mesh-1+i -max_side_of_mesh-1+j -max_side_of_mesh-1+k];
 tumour_GC.state_per_time(time)=0;
tumour_GC.num_of_cells_in_G1_per_time(time)=0;
 tumour_GC.num_of_cells_in_S_per_time(time)= 0;
tumour_GC.num_of_cells_in_G2_per_time(time)=0;
 tumour_GC.num_of_cells_in_M_per_time(time)=0;
 tumour_GC.num_of_cells_in_G0_per_time(time)=0;
tumour_GC.num_of_cells_in_N_per_time(time)=0;
tumour_GC.num_of_cells_in_necrotic_per_time(time)=0;
tumour_GC.num_of_all_cells_per_time(time)=0;
tumour_GC.dim_for_cells_radiated=0;
tumour_GC.cells_radiated=[0 0];      
tumour_GC.num_of_prolif_cells_per_time(time)=0;
tumour_GC.available_space_per_time(time)=0;
tumour_GC.prof_to_all(time)=0;
tumour_GC.rest_to_all(time)=0;
tumour_GC.dead_to_all(time)=0;
                    

tumour_GC.time_spent_in_G1_per_time=[];
tumour_GC.time_spent_in_S_per_time=[];
 tumour_GC.time_spent_in_G2_per_time=[];
 tumour_GC.time_spent_in_M_per_time=[];
tumour_GC.time_spent_in_G0_per_time=[];
 tumour_GC.time_spent_in_N_per_time=[];

 

   

