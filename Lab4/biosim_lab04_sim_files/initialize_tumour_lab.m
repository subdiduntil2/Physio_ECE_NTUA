function  tumour_in_GCs=initialize_tumour_lab(prolif_reg_exter_radii,rest_reg_exter_radii,necr_reg_exter_radii)

% It initializes the tumour
% It creates the tumour in space based on the radious of the 3 regions: prolferation, resting, necrotic from outside to inside.
% The tumour is space created as a struct: "tumour_in_GCs"
% Proliferating, Resting and Necrotic are two concentric rings and a sphere correspondingly

% prolif_reg_exter_radii cooresponds to the outer radii of the region of the proliferating GCs e.g. 8
% rest_reg_exter_radii cooresponds to the outer radii of the region of the resting GCs   e.g 5
% necr_reg_exter_radii cooresponds to the outer radii of the region of the necrotic GCs e.g 4


max_cell_density=10.^6; % we assume for all regions the maximum cell_density /mm.^3 (maximum number of cells in a GC) equal to 10.^6

init_cell_densit_in_prol_region=6*10.^5;
init_cell_densit_in_rest_region=6*10.^5;     %we give the initial cell_density in the 3 regions
init_cell_densit_in_necr_region=6*10.^5;

max_side_of_mesh=2*prolif_reg_exter_radii;  

% we assume maximum mesh side (for each of six directions equal to the double of prol_reg_exter_radii.
%for example: if prolif_reg_exter_radii=8, max_side_of_mesh=16. Then, the geometric mesh is the following x,y,z=-16..16
% and contains 33*33*33 geometric cells (GCs).                                             
                                           
                                         
T_C=6;       % the duration of the cell cycle, containing G1, S, G2 and M phases     
T_G1=2;
T_S=2;
T_G2=1;   
T_M=1;    %the duration of phases are assigned. Each time step in the programme corresponds to 6 hours, so e.g. the duration of G1 is 2*6=12 hours
T_G0=4;
T_N=8;

% The tumour is created as a struct with three dimensions
% When prolif_reg_exter_radii=8 then the tumour is created as a struct with 33*33*33 GCs, corresponding to GCs with coordinates x,y,z=-16..16

%Each GC has the following elements
%    coords .It containes the coordinates of the GC is space i.e tumour_in_GCs(1,1,1).coords=[-16 16 -16], or tumour_in_GCs(17,17,17).coords=[0 0 0]( when prolif_reg_exter_radii=8)
%    state_per_time .It indicates in the state of the GC: 0 corresponds to region out of tumour, 1 corresponds to the proliferating state, 2 corresponds to the resting state,3 corresponds to the necrotic state  i.e if tumour_in_GCs(17,17,17).state_per_time=[3 3 3 1] then the GC with coords 0 0 0 was in nevrotic in time 1 and in the proliferating in time 4
%    num_of_cells_in_G1_per_time   It shows the number of cells in G1 for the GC along time
%    num_of_cells_in_S_per_time      as above
%    num_of_cells_in_G2_per_time     as above
%    num_of_cells_in_M_per_time      as above
%    num_of_cells_in_G0_per_time     as above
%    num_of_cells_in_N_per_time      as above  (It shows the cells that are in N but are there due to the cycle and not due to radiation)
%    num_of_cells_in_necrotic_per_ti  as above  (It shows the cells that are in N due to the cycle and due to radiation)
%    num_of_all_cells_per_time        as above
%    dim_for_cells_radiated
%    cells_radiated
%    num_of_prolif_cells_per_time     as above (It shoes the cells in G1, S, G2, M)
%    available_space_per_time         corresponds to max_cell_density-num_of_all_cells
%    prof_to_all                      the percentage of the proligerating (G1,S,G2,M) to all
%    rest_to_all                      the percentage of the resting cells (G0) to all
%    dead_to_all                      the percentage of the dead (num_of_cells_in_necrotic) to all
%    time_spent_in_G1_per_time        it shows the time the cells in the GC have spent in G1 (it take values 1 to 2)
%    time_spent_in_S_per_time         as above
%    time_spent_in_G2_per_time        as above
%    time_spent_in_M_per_time         as above
%    time_spent_in_G0_per_time        as above
%    time_spent_in_N_per_time         as above 


%The following code initializes the tumour
for i=1:2*max_side_of_mesh+1;
    for j=1:2*max_side_of_mesh+1;
        for k=1:2*max_side_of_mesh+1;
            
            tumour_in_GCs(i,j,k).coords=[-max_side_of_mesh-1+i -max_side_of_mesh-1+j -max_side_of_mesh-1+k]; 
            dist_GC_center=sqrt(tumour_in_GCs(i,j,k).coords(1).^2+tumour_in_GCs(i,j,k).coords(2).^2+tumour_in_GCs(i,j,k).coords(3).^2);
            
            if (dist_GC_center<=necr_reg_exter_radii)
                
                tumour_in_GCs(i,j,k).state_per_time(1)=3; % 3 corresponds to necrotic region
                
                tumour_in_GCs(i,j,k).num_of_cells_in_G1_per_time(1)=T_G1/T_C*0.02*init_cell_densit_in_necr_region;
                tumour_in_GCs(i,j,k).num_of_cells_in_S_per_time(1)=T_S/T_C*0.02*init_cell_densit_in_necr_region;
                tumour_in_GCs(i,j,k).num_of_cells_in_G2_per_time(1)=T_G2/T_C*0.02*init_cell_densit_in_necr_region;  %we initialize the cell_numbers in various phases in the necrotic region 
                tumour_in_GCs(i,j,k).num_of_cells_in_M_per_time(1)=T_M/T_C*0.02*init_cell_densit_in_necr_region;
                tumour_in_GCs(i,j,k).num_of_cells_in_G0_per_time(1)=0.07*init_cell_densit_in_necr_region;                
                tumour_in_GCs(i,j,k).num_of_cells_in_N_per_time(1)=0.91*init_cell_densit_in_necr_region;     
                tumour_in_GCs(i,j,k).num_of_cells_in_necrotic_per_time(1)=tumour_in_GCs(i,j,k).num_of_cells_in_N_per_time(1);
                
                tumour_in_GCs(i,j,k).dim_for_cells_radiated=0;
                tumour_in_GCs(i,j,k).cells_radiated=[0 0];
                
                tumour_in_GCs(i,j,k).num_of_all_cells_per_time(1)=tumour_in_GCs(i,j,k).num_of_cells_in_G1_per_time(1)+tumour_in_GCs(i,j,k).num_of_cells_in_S_per_time(1)+tumour_in_GCs(i,j,k).num_of_cells_in_G2_per_time(1)+tumour_in_GCs(i,j,k).num_of_cells_in_G0_per_time(1)+tumour_in_GCs(i,j,k).num_of_cells_in_M_per_time(1)+ tumour_in_GCs(i,j,k).num_of_cells_in_N_per_time(1);
                tumour_in_GCs(i,j,k).num_of_prolif_cells_per_time(1)=tumour_in_GCs(i,j,k).num_of_cells_in_G1_per_time(1)+tumour_in_GCs(i,j,k).num_of_cells_in_S_per_time(1)+tumour_in_GCs(i,j,k).num_of_cells_in_G2_per_time(1)+tumour_in_GCs(i,j,k).num_of_cells_in_M_per_time(1);
                tumour_in_GCs(i,j,k).available_space_per_time(1)=max(max_cell_density-tumour_in_GCs(i,j,k).num_of_all_cells_per_time(1),0);
                
                               
                tumour_in_GCs(i,j,k).prof_to_all(1)= tumour_in_GCs(i,j,k).num_of_prolif_cells_per_time(1)/tumour_in_GCs(i,j,k).num_of_all_cells_per_time(1);
                tumour_in_GCs(i,j,k).rest_to_all(1)= tumour_in_GCs(i,j,k).num_of_cells_in_G0_per_time(1)/tumour_in_GCs(i,j,k).num_of_all_cells_per_time(1);
                tumour_in_GCs(i,j,k).dead_to_all(1)= tumour_in_GCs(i,j,k).num_of_cells_in_N_per_time(1)/tumour_in_GCs(i,j,k).num_of_all_cells_per_time(1);
               
                tumour_in_GCs(i,j,k).time_spent_in_G1_per_time(1)=ceil(rand*T_G1);
                tumour_in_GCs(i,j,k).time_spent_in_S_per_time(1)=ceil(rand*T_S);
                tumour_in_GCs(i,j,k).time_spent_in_G2_per_time(1)=ceil(rand*T_G2);  %when it is e.g. 3 it means that in the next run it is the cell's 4th hour.
                tumour_in_GCs(i,j,k).time_spent_in_M_per_time(1)=ceil(rand*T_M);
                tumour_in_GCs(i,j,k).time_spent_in_G0_per_time(1)=ceil(rand*T_G0);
                tumour_in_GCs(i,j,k).time_spent_in_N_per_time(1)=ceil(rand*T_N);
                                
            elseif (dist_GC_center<=rest_reg_exter_radii)  
            
                tumour_in_GCs(i,j,k).state_per_time(1)=2; % 2 corresponds to resting region
                
                tumour_in_GCs(i,j,k).num_of_cells_in_G1_per_time(1)=T_G1/T_C*0.1*init_cell_densit_in_rest_region;
                tumour_in_GCs(i,j,k).num_of_cells_in_S_per_time(1)=T_S/T_C*0.1*init_cell_densit_in_rest_region;
                tumour_in_GCs(i,j,k).num_of_cells_in_G2_per_time(1)=T_G2/T_C*0.1*init_cell_densit_in_rest_region;  %we initialize the cell_numbers in various phases in the resting region 
                tumour_in_GCs(i,j,k).num_of_cells_in_M_per_time(1)=T_M/T_C*0.1*init_cell_densit_in_rest_region;
                tumour_in_GCs(i,j,k).num_of_cells_in_G0_per_time(1)=0.8*init_cell_densit_in_rest_region;                
                tumour_in_GCs(i,j,k).num_of_cells_in_N_per_time(1)=0.1*init_cell_densit_in_rest_region;
                tumour_in_GCs(i,j,k).num_of_cells_in_necrotic_per_time(1)=tumour_in_GCs(i,j,k).num_of_cells_in_N_per_time(1);
                
                tumour_in_GCs(i,j,k).dim_for_cells_radiated=0;
                tumour_in_GCs(i,j,k).cells_radiated=[0 0];
                
                tumour_in_GCs(i,j,k).num_of_all_cells_per_time(1)=tumour_in_GCs(i,j,k).num_of_cells_in_G1_per_time(1)+tumour_in_GCs(i,j,k).num_of_cells_in_S_per_time(1)+tumour_in_GCs(i,j,k).num_of_cells_in_G2_per_time(1)+tumour_in_GCs(i,j,k).num_of_cells_in_G0_per_time(1)+tumour_in_GCs(i,j,k).num_of_cells_in_M_per_time(1)+ tumour_in_GCs(i,j,k).num_of_cells_in_N_per_time(1);
                tumour_in_GCs(i,j,k).num_of_prolif_cells_per_time(1)=tumour_in_GCs(i,j,k).num_of_cells_in_G1_per_time(1)+tumour_in_GCs(i,j,k).num_of_cells_in_S_per_time(1)+tumour_in_GCs(i,j,k).num_of_cells_in_G2_per_time(1)+tumour_in_GCs(i,j,k).num_of_cells_in_M_per_time(1);
                tumour_in_GCs(i,j,k).available_space_per_time(1)=max(max_cell_density-tumour_in_GCs(i,j,k).num_of_all_cells_per_time(1),0);
                
                tumour_in_GCs(i,j,k).prof_to_all(1)= tumour_in_GCs(i,j,k).num_of_prolif_cells_per_time(1)/tumour_in_GCs(i,j,k).num_of_all_cells_per_time(1);
                tumour_in_GCs(i,j,k).rest_to_all(1)= tumour_in_GCs(i,j,k).num_of_cells_in_G0_per_time(1)/tumour_in_GCs(i,j,k).num_of_all_cells_per_time(1);
                tumour_in_GCs(i,j,k).dead_to_all(1)= tumour_in_GCs(i,j,k).num_of_cells_in_N_per_time(1)/tumour_in_GCs(i,j,k).num_of_all_cells_per_time(1);
               
                
                tumour_in_GCs(i,j,k).time_spent_in_G1_per_time(1)=ceil(rand*T_G1);
                tumour_in_GCs(i,j,k).time_spent_in_S_per_time(1)=ceil(rand*T_S);
                tumour_in_GCs(i,j,k).time_spent_in_G2_per_time(1)=ceil(rand*T_G2);  %when it is e.g. 3 it means that in the next run it is the cell's 4th hour.
                tumour_in_GCs(i,j,k).time_spent_in_M_per_time(1)=ceil(rand*T_M);
                tumour_in_GCs(i,j,k).time_spent_in_G0_per_time(1)=ceil(rand*T_G0);
                tumour_in_GCs(i,j,k).time_spent_in_N_per_time(1)=ceil(rand*T_N);
                
            elseif (dist_GC_center<=prolif_reg_exter_radii)  
            
                tumour_in_GCs(i,j,k).state_per_time(1)=1; % 1 corresponds to proliferating region
                
                tumour_in_GCs(i,j,k).num_of_cells_in_G1_per_time(1)=T_G1/T_C*0.8*init_cell_densit_in_prol_region;
                tumour_in_GCs(i,j,k).num_of_cells_in_S_per_time(1)=T_S/T_C*0.8*init_cell_densit_in_prol_region;
                tumour_in_GCs(i,j,k).num_of_cells_in_G2_per_time(1)=T_G2/T_C*0.8*init_cell_densit_in_prol_region;  %we initialize the cell_numbers in various phases in the proliferating region 
                tumour_in_GCs(i,j,k).num_of_cells_in_M_per_time(1)=T_M/T_C*0.8*init_cell_densit_in_prol_region;
                tumour_in_GCs(i,j,k).num_of_cells_in_G0_per_time(1)=0.18*init_cell_densit_in_prol_region;                
                tumour_in_GCs(i,j,k).num_of_cells_in_N_per_time(1)=0.02*init_cell_densit_in_prol_region;
                tumour_in_GCs(i,j,k).num_of_cells_in_necrotic_per_time(1)=tumour_in_GCs(i,j,k).num_of_cells_in_N_per_time(1);
                
                tumour_in_GCs(i,j,k).dim_for_cells_radiated=0;
                tumour_in_GCs(i,j,k).cells_radiated=[0 0];
                
                tumour_in_GCs(i,j,k).num_of_all_cells_per_time(1)=tumour_in_GCs(i,j,k).num_of_cells_in_G1_per_time(1)+tumour_in_GCs(i,j,k).num_of_cells_in_S_per_time(1)+tumour_in_GCs(i,j,k).num_of_cells_in_G2_per_time(1)+tumour_in_GCs(i,j,k).num_of_cells_in_G0_per_time(1)+tumour_in_GCs(i,j,k).num_of_cells_in_M_per_time(1)+ tumour_in_GCs(i,j,k).num_of_cells_in_N_per_time(1);
                tumour_in_GCs(i,j,k).num_of_prolif_cells_per_time(1)=tumour_in_GCs(i,j,k).num_of_cells_in_G1_per_time(1)+tumour_in_GCs(i,j,k).num_of_cells_in_S_per_time(1)+tumour_in_GCs(i,j,k).num_of_cells_in_G2_per_time(1)+tumour_in_GCs(i,j,k).num_of_cells_in_M_per_time(1);
                tumour_in_GCs(i,j,k).available_space_per_time(1)=max(max_cell_density-tumour_in_GCs(i,j,k).num_of_all_cells_per_time(1),0);
                
                tumour_in_GCs(i,j,k).prof_to_all(1)= tumour_in_GCs(i,j,k).num_of_prolif_cells_per_time(1)/tumour_in_GCs(i,j,k).num_of_all_cells_per_time(1);
                tumour_in_GCs(i,j,k).rest_to_all(1)= tumour_in_GCs(i,j,k).num_of_cells_in_G0_per_time(1)/tumour_in_GCs(i,j,k).num_of_all_cells_per_time(1);
                tumour_in_GCs(i,j,k).dead_to_all(1)= tumour_in_GCs(i,j,k).num_of_cells_in_N_per_time(1)/tumour_in_GCs(i,j,k).num_of_all_cells_per_time(1);
               
                tumour_in_GCs(i,j,k).time_spent_in_G1_per_time(1)=ceil(rand*T_G1);
                tumour_in_GCs(i,j,k).time_spent_in_S_per_time(1)=ceil(rand*T_S);
                tumour_in_GCs(i,j,k).time_spent_in_G2_per_time(1)=ceil(rand*T_G2);  %when it is e.g. 3 it means that in the next run it is the cell's 4th hour.
                tumour_in_GCs(i,j,k).time_spent_in_M_per_time(1)=ceil(rand*T_M);
                tumour_in_GCs(i,j,k).time_spent_in_G0_per_time(1)=ceil(rand*T_G0);
                tumour_in_GCs(i,j,k).time_spent_in_N_per_time(1)=ceil(rand*T_N);
                
            else    
                
                tumour_in_GCs(i,j,k).state_per_time(1)=0; % 0 corresponds to out of tumour region
                                
                tumour_in_GCs(i,j,k).num_of_cells_in_G1_per_time(1)=0;
                tumour_in_GCs(i,j,k).num_of_cells_in_S_per_time(1)=0;
                tumour_in_GCs(i,j,k).num_of_cells_in_G2_per_time(1)=0;  %we initialize the cell_numbers in various phases in the proliferating region 
                tumour_in_GCs(i,j,k).num_of_cells_in_M_per_time(1)=0;
                tumour_in_GCs(i,j,k).num_of_cells_in_G0_per_time(1)=0;                
                tumour_in_GCs(i,j,k).num_of_cells_in_N_per_time(1)=0;                   
                tumour_in_GCs(i,j,k).num_of_cells_in_necrotic_per_time(1)=0;
                tumour_in_GCs(i,j,k).num_of_all_cells_per_time(1)=0;    
                
                tumour_in_GCs(i,j,k).dim_for_cells_radiated=0;
                tumour_in_GCs(i,j,k).cells_radiated=[0 0];
                
                
                tumour_in_GCs(i,j,k).num_of_prolif_cells_per_time(1)=0;
                tumour_in_GCs(i,j,k).available_space_per_time(1)=max_cell_density;
                
                tumour_in_GCs(i,j,k).prof_to_all(1)= 0;
                tumour_in_GCs(i,j,k).rest_to_all(1)= 0;
                tumour_in_GCs(i,j,k).dead_to_all(1)= 0;
            end
            
        end
    end
end


% After the above code has run we can access for example the GC in the centre by typing: tumour_in_GCs(17,17,17)
    
    