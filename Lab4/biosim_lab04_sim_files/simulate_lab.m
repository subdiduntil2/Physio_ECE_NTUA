%This function simulates a tumour's growth and its response to
%radiotherapy..

%Before running the current function the following variables must be in the
%workspace:
%  a variable tumour_in_GCs that is a struct and contains the
%                    initialization of the tumour. you can set this
%                    variable by typing
%                    tumour_in_GCs=initialize_tumour_lab(8,5,4)
%  two variables :a","b" that contain the radiosensitivity parameters of the
%                    tumour. You can set these variables by typing for example a=[0.6 0.5 0.4] and b=[0.06 0.05 0.04]
% a variable "radiotherapy_scheme" that containes the radiotherapy_scheme.
%                   You can set this variable by typing radiotherapy_scheme=[0 2 0 0 0 2...]
%                   or inserting the corresping mat files provided to you
% a variable duration that is the duration of the simulation in time steps.You can type
% e.g. duration=57 (for simulation for two weeks!!) or [anything, duration]=size(radiotherapy_scheme)



[m,m,m]=size(tumour_in_GCs); % The length of "tumour_in_GCs" is assigned to m

max_cell_density=10.^6; %The max_cell_density is assigned to the corresponding variable

T_C=6;      
T_G1=2;
T_S=2;
T_G2=1;   %The durations in time steps of the G1,S,G2,M,G0,N phases are assigned
T_M=1;   
T_G0=4;
T_N=8;


%a=[a_P a_S a_Go]
a_P=a(1);
a_S=a(2);       % we extract from vector a the values for the parameters a_P,a_S,a_G0 corresponding to sensitivity parameter in (G1,G2,M), S and G0 
a_G0=a(3);


%b=[b_P b_S b_Go]
b_P=b(1);
b_S=b(2);       % we extract from vector b the values for the parameters b_P,b_S,b_G0 corresponding to to sensitivity parameter in (G1,G2,M), S and G0 
b_G0=b(3);


%we assign the propabilities of transiting from M to G1,G0 and from G0 to G1,N for all states (the first element is for a GC in proliferating state, the second element is for a GC in resting state, the third element is for a GC in necrotic state)  
Probs_M_to_G1=[0.7 0.3 0.2];  % Probs_M_to_G1(1) is for proliferationg state, Probs_M_to_G1(2) is for resting state, Probs_M_to_G1(3) is for necrotic state
Probs_M_to_G0=[0.3 0.7 0.8];
Probs_G0_to_G1=[0.8 0.1 0.02];
Probs_G0_to_N=[0.2 0.9 0.98];


for time=2:duration  % time scanning begins at the step time=2 because e.g tumour_in_GCs(i,j,k).num_of_cells_in_G1_per_time(1) is for the initialized state
       
   fprintf('Time step= %d\n', time)  % the time step is displayed so as we can see the progress of the simulation procedure
 
        D=radiotherapy_scheme(time);   % D is the radiation at the current time step, D can be zero or not
        for i=1:m       % 1st space scanning, All GCs are scanned and the changes in the number of cells due to radiation are calculaeted       
            for j=1:m
                for k=1:m
                                                  
                     tumour_in_GCs(i,j,k).num_of_cells_in_G1_per_time(time)=tumour_in_GCs(i,j,k).num_of_cells_in_G1_per_time(time-1)*exp(-(a_P*D+b_P*D.^2));
                     tumour_in_GCs(i,j,k).num_of_cells_in_S_per_time(time)=tumour_in_GCs(i,j,k).num_of_cells_in_S_per_time(time-1)*exp(-(a_S*D+b_S*D.^2));  % ADD YOUR CODE AS ABOVE
                     tumour_in_GCs(i,j,k).num_of_cells_in_G2_per_time(time)=tumour_in_GCs(i,j,k).num_of_cells_in_G2_per_time(time-1)*exp(-(a_P*D+b_P*D.^2));  % ADD YOUR CODE AS ABOVE                                                                              % in case of radiation some of the cells die in phases G1,G2,S,M,G0. if D=0 no cells die.
                     tumour_in_GCs(i,j,k).num_of_cells_in_M_per_time(time)=tumour_in_GCs(i,j,k).num_of_cells_in_M_per_time(time-1)*exp(-(a_P*D+b_P*D.^2));  % ADD YOUR CODE AS ABOVE                                                                               % use the values of a_P,a_S,a_G0, a_P,a_S,a_G0 as defined in the beginning of this file
                     tumour_in_GCs(i,j,k).num_of_cells_in_G0_per_time(time)=tumour_in_GCs(i,j,k).num_of_cells_in_G0_per_time(time-1)*exp(-(a_G0*D+b_G0*D.^2)); % ADD YOUR CODE AS ABOVE
                     
                     if D~=0
                         tumour_in_GCs(i,j,k).dim_for_cells_radiated=tumour_in_GCs(i,j,k).dim_for_cells_radiated+1;    %if D~=0, then the percentage of dying cells due to radiation will disappear after 12 hours(=2 * duration of the cell cycle containing G1,S, G2,M)
                         tumour_in_GCs(i,j,k).cells_radiated(tumour_in_GCs(i,j,k).dim_for_cells_radiated,1:2)=[tumour_in_GCs(i,j,k).num_of_cells_in_G1_per_time(time-1)-tumour_in_GCs(i,j,k).num_of_cells_in_G1_per_time(time)+tumour_in_GCs(i,j,k).num_of_cells_in_G2_per_time(time-1)-tumour_in_GCs(i,j,k).num_of_cells_in_G2_per_time(time)+tumour_in_GCs(i,j,k).num_of_cells_in_S_per_time(time-1)-tumour_in_GCs(i,j,k).num_of_cells_in_S_per_time(time)+tumour_in_GCs(i,j,k).num_of_cells_in_M_per_time(time-1)-tumour_in_GCs(i,j,k).num_of_cells_in_M_per_time(time)+tumour_in_GCs(i,j,k).num_of_cells_in_G0_per_time(time-1)-tumour_in_GCs(i,j,k).num_of_cells_in_G0_per_time(time) time+13];
                     end
                     
                     tumour_in_GCs(i,j,k).num_of_cells_in_N_per_time(time)=tumour_in_GCs(i,j,k).num_of_cells_in_N_per_time(time-1);  % the cell in necrosis stay the same as the previous time step for now
                     
                     tumour_in_GCs(i,j,k).num_of_all_cells_per_time(time)=tumour_in_GCs(i,j,k).num_of_cells_in_G1_per_time(time)+tumour_in_GCs(i,j,k).num_of_cells_in_S_per_time(time)+tumour_in_GCs(i,j,k).num_of_cells_in_G2_per_time(time)+tumour_in_GCs(i,j,k).num_of_cells_in_G0_per_time(time)+tumour_in_GCs(i,j,k).num_of_cells_in_M_per_time(time)+ tumour_in_GCs(i,j,k).num_of_cells_in_N_per_time(time);
                     for radiat_cases=1:tumour_in_GCs(i,j,k).dim_for_cells_radiated
                     tumour_in_GCs(i,j,k).num_of_all_cells_per_time(time)=tumour_in_GCs(i,j,k).num_of_all_cells_per_time(time)+tumour_in_GCs(i,j,k).cells_radiated(radiat_cases,1);   % all the cells of the GC (G1,G2,S,M, GO, N and the dead due to radiation) are calculated
                     end    
                     
                 end  % for k=..
            end   %for j=1...   End of 1st space scanning (changes due to radiation)
        end   %for i=1...   
            

            
   
     for i=1:m          %2nd space scanning, All GCs are scanned and the phases transitions are made..Also (in the very beginning) the loss of cells in G1,S, G2, M due to apoptosis is calculated  
            for j=1:m
                for k=1:m
                   
                    
                    
                  if tumour_in_GCs(i,j,k).num_of_all_cells_per_time(time)~=0
                                      
                    
                    %loss of cells due to apoptosis  
                    tumour_in_GCs(i,j,k).num_of_cells_in_G1_per_time(time)=tumour_in_GCs(i,j,k).num_of_cells_in_G1_per_time(time)*0.999999; % 1/1000000 are lost due to apoptosis
                    tumour_in_GCs(i,j,k).num_of_cells_in_S_per_time(time)=tumour_in_GCs(i,j,k).num_of_cells_in_S_per_time(time)*0.999999;
                    tumour_in_GCs(i,j,k).num_of_cells_in_G2_per_time(time)=tumour_in_GCs(i,j,k).num_of_cells_in_G2_per_time(time)*0.999999;
                    tumour_in_GCs(i,j,k).num_of_cells_in_M_per_time(time)=tumour_in_GCs(i,j,k).num_of_cells_in_M_per_time(time)*0.999999;                 %ADD YOUR CODE AS ABOVE
                    tumour_in_GCs(i,j,k).num_of_cells_in_G0_per_time(time)=tumour_in_GCs(i,j,k).num_of_cells_in_G0_per_time(time)*0.999999;
                       
                    
                    [tumour_in_GCs(i,j,k).dim_for_cells_radiated, ka]=size(tumour_in_GCs(i,j,k).cells_radiated);
                    for radiat_cases=1:tumour_in_GCs(i,j,k).dim_for_cells_radiated
                           if tumour_in_GCs(i,j,k).cells_radiated(radiat_cases,2)==time
                           tumour_in_GCs(i,j,k).cells_radiated(radiat_cases,:)=[0 0];       % if 12 hours have passed from radiation for one of the radiation cases occured until now the corresponding dead cells disappear
                           end
                    end   
                    
                    
                    %calculation of cells for from phase to phase transition   
                    if  tumour_in_GCs(i,j,k).time_spent_in_G1_per_time(time-1)==T_G1
                        cells_from_G1_to_S=tumour_in_GCs(i,j,k).num_of_cells_in_G1_per_time(time);    % cells that transit from G1 to S
                                     
                    else
                        cells_from_G1_to_S=0;
                    end

                    if  tumour_in_GCs(i,j,k).time_spent_in_S_per_time(time-1)==T_S
                        cells_from_S_to_G2=tumour_in_GCs(i,j,k).num_of_cells_in_S_per_time(time);    %cells that transit from S to G2
                                     
                    else
                        cells_from_S_to_G2=0;
                    end
                        
                    if  tumour_in_GCs(i,j,k).time_spent_in_G2_per_time(time-1)==T_G2
                        cells_from_G2_to_M=tumour_in_GCs(i,j,k).num_of_cells_in_G2_per_time(time);    %cells that transit from G2 to M
                                     
                    else
                        cells_from_G2_to_M=0;
                    end
                       
                     if  tumour_in_GCs(i,j,k).time_spent_in_M_per_time(time-1)==T_M
                        tumour_in_GCs(i,j,k).num_of_cells_in_M_per_time(time)=2*tumour_in_GCs(i,j,k).num_of_cells_in_M_per_time(time);
                        cells_from_M_to_G1=Probs_M_to_G1(tumour_in_GCs(i,j,k).state_per_time(time-1))*tumour_in_GCs(i,j,k).num_of_cells_in_M_per_time(time);    %cells that transit from M to G1 or G0        %if M phase has ended, firstly the cells in M get doubled and then the transitions are made                     %USE THE PROPABILITIES FOR THE TRANSITIONS FROM PHASE TO PHASE WHERE APPLICABLE (these propabilities have been defined at the start of this file: Probs_M_to_G1 etc.)
                        cells_from_M_to_G0=Probs_M_to_G0(tumour_in_GCs(i,j,k).state_per_time(time-1))*tumour_in_GCs(i,j,k).num_of_cells_in_M_per_time(time);                                                                                                                                                                                                                                            %the propability depends on state the GC was in the previous time step. Get for example the correct propability by typing: Probs_M_to_G1(tumour_in_GCs(i,j,k).state_per_time(time-1))
                                     
                    else
                        cells_from_M_to_G1=0;
                        cells_from_M_to_G0=0;
                    end
                    
                    if  tumour_in_GCs(i,j,k).time_spent_in_G0_per_time(time-1)==T_G0
                        cells_from_G0_to_G1=Probs_G0_to_G1(tumour_in_GCs(i,j,k).state_per_time(time-1))*tumour_in_GCs(i,j,k).num_of_cells_in_G0_per_time(time);    %cells that transit from G0 to G1 or N       %if M phase has ended, firstly the cells in M get doubled and then the transitions are made                     %USE THE PROPABILITIES FOR THE TRANSITIONS FROM PHASE TO PHASE WHERE APPLICABLE (these propabilities have been defined at the start of this file: Probs_M_to_G1 etc.)
                        cells_from_G0_to_N=Probs_G0_to_N(tumour_in_GCs(i,j,k).state_per_time(time-1))*tumour_in_GCs(i,j,k).num_of_cells_in_G0_per_time(time);                                                                                                                                                                                                                                            %the propability depends on state the GC was in the previous time step. Get for example the correct propability by typing: Probs_M_to_G1(tumour_in_GCs(i,j,k).state_per_time(time-1))
                                     
                    else
                        cells_from_G0_to_G1=0;
                        cells_from_G0_to_N=0;
                    end

                    %calculation of cells for from phase to phase transition   
                    if  tumour_in_GCs(i,j,k).time_spent_in_N_per_time(time-1)==T_N
                        cells_from_N_to_lysis=tumour_in_GCs(i,j,k).num_of_cells_in_N_per_time(time);    % cells that transit from G1 to S
                                     
                    else
                        cells_from_N_to_lysis=0;
                    end
                    
                    % ADD YOUR CODE AS ABOVE                                                            %cells that transit from N to lysis and get lost
                        
                                       
                    tumour_in_GCs(i,j,k).num_of_cells_in_G1_per_time(time)=tumour_in_GCs(i,j,k).num_of_cells_in_G1_per_time(time)-cells_from_G1_to_S+cells_from_M_to_G1+cells_from_G0_to_G1;
                    tumour_in_GCs(i,j,k).num_of_cells_in_S_per_time(time)= tumour_in_GCs(i,j,k).num_of_cells_in_S_per_time(time)-cells_from_S_to_G2+cells_from_G1_to_S;   % ADD YOUR CODE AS ABOVE
                    tumour_in_GCs(i,j,k).num_of_cells_in_G2_per_time(time)=tumour_in_GCs(i,j,k).num_of_cells_in_G2_per_time(time)-cells_from_G2_to_M+cells_from_S_to_G2;  % ADD YOUR CODE AS ABOVE
                    tumour_in_GCs(i,j,k).num_of_cells_in_M_per_time(time)= tumour_in_GCs(i,j,k).num_of_cells_in_M_per_time(time) - cells_from_M_to_G0 -cells_from_M_to_G1 +cells_from_G2_to_M;   % ADD YOUR CODE AS ABOVE         %BASED ON THE CELL CYCLE SCHEME          %this code find the cells in the GC left after apoptosis and phase to phase transition (the changes to due to radiation have already been calculated in the first space scanning) 
                    tumour_in_GCs(i,j,k).num_of_cells_in_G0_per_time(time)=tumour_in_GCs(i,j,k).num_of_cells_in_G0_per_time(time) -cells_from_G0_to_N -cells_from_G0_to_G1 + cells_from_M_to_G0;    % ADD YOUR CODE AS ABOVE
                    tumour_in_GCs(i,j,k).num_of_cells_in_N_per_time(time)=tumour_in_GCs(i,j,k).num_of_cells_in_N_per_time(time) - cells_from_N_to_lysis + cells_from_G0_to_N;   % ADD YOUR CODE AS ABOVE                 
                    
                    tumour_in_GCs(i,j,k).num_of_all_cells_per_time(time)=tumour_in_GCs(i,j,k).num_of_cells_in_G1_per_time(time)+tumour_in_GCs(i,j,k).num_of_cells_in_S_per_time(time)+tumour_in_GCs(i,j,k).num_of_cells_in_G2_per_time(time)+tumour_in_GCs(i,j,k).num_of_cells_in_M_per_time(time)+tumour_in_GCs(i,j,k).num_of_cells_in_G0_per_time(time)+tumour_in_GCs(i,j,k).num_of_cells_in_N_per_time(time);
                    for radiat_cases=1:tumour_in_GCs(i,j,k).dim_for_cells_radiated
                     tumour_in_GCs(i,j,k).num_of_all_cells_per_time(time)=tumour_in_GCs(i,j,k).num_of_all_cells_per_time(time)+tumour_in_GCs(i,j,k).cells_radiated(radiat_cases,1);  %in the last 3 lines:  all final cells are calculated for a GC including: cells in G1, S, G2, M, G0, N and cells dead to radiation
                    end
                    
                    tumour_in_GCs(i,j,k).num_of_prolif_cells_per_time(time)=tumour_in_GCs(i,j,k).num_of_cells_in_G1_per_time(time)+tumour_in_GCs(i,j,k).num_of_cells_in_S_per_time(time)+tumour_in_GCs(i,j,k).num_of_cells_in_G2_per_time(time)+tumour_in_GCs(i,j,k).num_of_cells_in_M_per_time(time);  %for each GC all cells that are proligerating (G1, S, G2, M ) are calculated  
                    tumour_in_GCs(i,j,k).available_space_per_time(time)=max(max_cell_density-tumour_in_GCs(i,j,k).num_of_all_cells_per_time(time),0); % the available space is calculated for each gc
                    
                    %time in each phase for each geometric cell becomes 1 greater if the phase duration has not ended or becomes 1 if it starts now
                    if  tumour_in_GCs(i,j,k).time_spent_in_G1_per_time(time-1)==T_G1
                        tumour_in_GCs(i,j,k).time_spent_in_G1_per_time(time)=1;
                    else   
                        tumour_in_GCs(i,j,k).time_spent_in_G1_per_time(time)=tumour_in_GCs(i,j,k).time_spent_in_G1_per_time(time-1)+1;
                    end
                    
                    if  tumour_in_GCs(i,j,k).time_spent_in_S_per_time(time-1)==T_S
                        tumour_in_GCs(i,j,k).time_spent_in_S_per_time(time)=1;
                    else   
                        tumour_in_GCs(i,j,k).time_spent_in_S_per_time(time)=tumour_in_GCs(i,j,k).time_spent_in_S_per_time(time-1)+1;
                    end
                    
                    if  tumour_in_GCs(i,j,k).time_spent_in_G2_per_time(time-1)==T_G2
                        tumour_in_GCs(i,j,k).time_spent_in_G2_per_time(time)=1;
                    else   
                        tumour_in_GCs(i,j,k).time_spent_in_G2_per_time(time)=tumour_in_GCs(i,j,k).time_spent_in_G2_per_time(time-1)+1;
                    end
                    
                    if  tumour_in_GCs(i,j,k).time_spent_in_M_per_time(time-1)==T_M
                        tumour_in_GCs(i,j,k).time_spent_in_M_per_time(time)=1;
                    else   
                        tumour_in_GCs(i,j,k).time_spent_in_M_per_time(time)=tumour_in_GCs(i,j,k).time_spent_in_M_per_time(time-1)+1;
                    end
                    
                    if  tumour_in_GCs(i,j,k).time_spent_in_G0_per_time(time-1)==T_G0
                        tumour_in_GCs(i,j,k).time_spent_in_G0_per_time(time)=1;
                    else   
                        tumour_in_GCs(i,j,k).time_spent_in_G0_per_time(time)=tumour_in_GCs(i,j,k).time_spent_in_G0_per_time(time-1)+1;
                    end
                    
                    if  tumour_in_GCs(i,j,k).time_spent_in_N_per_time(time-1)==T_N
                        tumour_in_GCs(i,j,k).time_spent_in_N_per_time(time)=1;
                    else   
                        tumour_in_GCs(i,j,k).time_spent_in_N_per_time(time)=tumour_in_GCs(i,j,k).time_spent_in_N_per_time(time-1)+1;
                    end
                    
                 
                 else
                    tumour_in_GCs(i,j,k).num_of_all_cells_per_time(time)=0;        
                    tumour_in_GCs(i,j,k).num_of_prolif_cells_per_time(time)=0;    
                    tumour_in_GCs(i,j,k).available_space_per_time(time)=max_cell_density;
                                                                                                   % assignments for GCs that are empty
                    tumour_in_GCs(i,j,k).prof_to_all(time)= 0;
                    tumour_in_GCs(i,j,k).rest_to_all(time)= 0;
                    tumour_in_GCs(i,j,k).dead_to_all(time)= 0;    
                    
                 end  % for if ( tumour_in_GCs(i,j,k).num_of_all_cells_per_time(time)~=0)
                 
                             
                              
                             
                end  % for i=1..
            end    % for j=1..              End of 2nd space scanning  (phase to phase transition and  apoptosis)
        end   % for k=1..

           
        for i=1:m     % 3rd space scanning .It scanns each GC and when it containes cells less than 0.5*max_cell_density, the cells of the GC are transferred to neighbor GCs and that GC disappears
           for j=1:m
              for k=1:m
                
                  if (tumour_in_GCs(i,j,k).num_of_all_cells_per_time(time)<=0.5*max_cell_density) &(tumour_in_GCs(i,j,k).num_of_all_cells_per_time(time)~=0)
                      
                     
                      neighbourhood_availability=zeros(1,1);
                  
                      
                      num_of_neighbors=0;
                      
                      neighbors_i_s=zeros(1,1);
                      neighbors_j_s=zeros(1,1);
                      neighbors_k_s=zeros(1,1);
                      
                      % with the 6 following if..end command blocks the neighbors in the 6  direction (+-x, +-y, +-z) are found and if available for taking some content their availability is calculated  
                      if (i-1>=1) & (tumour_in_GCs(i-1,j,k).num_of_all_cells_per_time(time)~=0) & (tumour_in_GCs(i-1,j,k).available_space_per_time(time)>0)  % the neighbor in the -x direction
                          num_of_neighbors=num_of_neighbors+1;
                          neighbors_i_s(num_of_neighbors)=-1;
                          neighbors_j_s(num_of_neighbors)=0;
                          neighbors_k_s(num_of_neighbors)=0;
                          neighbourhood_availability(num_of_neighbors)=tumour_in_GCs(i-1,j,k).available_space_per_time(time);   % availability of the neighbor is the available space of the neighbor 
                          
                      end
     
                      if (j-1>=1) &(tumour_in_GCs(i,j-1,k).num_of_all_cells_per_time(time)~=0) & (tumour_in_GCs(i,j-1,k).available_space_per_time(time)>0)  % the neighbor in the -y direction
                          num_of_neighbors=num_of_neighbors+1;
                          neighbors_i_s(num_of_neighbors)=0;
                          neighbors_j_s(num_of_neighbors)=-1;
                          neighbors_k_s(num_of_neighbors)=0;
                          neighbourhood_availability(num_of_neighbors)=tumour_in_GCs(i,j-1,k).available_space_per_time(time);   
                      end
         
                      if (k-1>=1) & (tumour_in_GCs(i,j,k-1).num_of_all_cells_per_time(time)~=0) & (tumour_in_GCs(i,j,k-1).available_space_per_time(time)>0)  % the neighbor in the -z direction
                          num_of_neighbors=num_of_neighbors+1;
                          neighbors_i_s(num_of_neighbors)=0;
                          neighbors_j_s(num_of_neighbors)=0;
                          neighbors_k_s(num_of_neighbors)=-1;
                          neighbourhood_availability(num_of_neighbors)=tumour_in_GCs(i,j,k-1).available_space_per_time(time);     
                      end
        
                      if (i+1<=m) & (tumour_in_GCs(i+1,j,k).num_of_all_cells_per_time(time)~=0) & (tumour_in_GCs(i+1,j,k).available_space_per_time(time)>0)  % the neighbor in the +x direction
                          num_of_neighbors=num_of_neighbors+1;
                          neighbors_i_s(num_of_neighbors)=1;
                          neighbors_j_s(num_of_neighbors)=0;
                          neighbors_k_s(num_of_neighbors)=0;
                          neighbourhood_availability(num_of_neighbors)=tumour_in_GCs(i+1,j,k).available_space_per_time(time);    
                      end
        
                      if (j+1<=m) &(tumour_in_GCs(i,j+1,k).num_of_all_cells_per_time(time)~=0) & (tumour_in_GCs(i,j+1,k).available_space_per_time(time)>0) % the neighbor in the +y direction
                          num_of_neighbors=num_of_neighbors+1;
                          neighbors_i_s(num_of_neighbors)=0;
                          neighbors_j_s(num_of_neighbors)=1;
                          neighbors_k_s(num_of_neighbors)=0;
                          neighbourhood_availability(num_of_neighbors)=tumour_in_GCs(i,j+1,k).available_space_per_time(time);   
                      end
        
                      if (k+1<=m) &(tumour_in_GCs(i,j,k+1).num_of_all_cells_per_time(time)~=0) & (tumour_in_GCs(i,j,k+1).available_space_per_time(time)>0)   % the neighbor in the +z direction
                          num_of_neighbors=num_of_neighbors+1;
                          neighbors_i_s(num_of_neighbors)=0;
                          neighbors_j_s(num_of_neighbors)=0;
                          neighbors_k_s(num_of_neighbors)=1;
                          neighbourhood_availability(num_of_neighbors)=tumour_in_GCs(i,j,k+1).available_space_per_time(time);  
                      end
                      
                      for neighbor=1:num_of_neighbors  % this for loop transfers the content of the GC to the suitable neighbors depending on their available space, and changes the times_spent in the phases depending on the cells exchange
                      
                          fraction=neighbourhood_availability(neighbor)/sum(neighbourhood_availability);
                          
                          %the following block trasfers the cells in G1 phase from the GC to the neighbor and changed finds the suitable time_spent_in_G1 for the neighbor 
                          divider=(fraction*tumour_in_GCs(i,j,k).num_of_cells_in_G1_per_time(time)+tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).num_of_cells_in_G1_per_time(time));
                          if (divider~=0)   
                          tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).time_spent_in_G1_per_time(time)=round((tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).num_of_cells_in_G1_per_time(time)*tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).time_spent_in_G1_per_time(time)+fraction*tumour_in_GCs(i,j,k).num_of_cells_in_G1_per_time(time)*tumour_in_GCs(i,j,k).time_spent_in_G1_per_time(time))/divider);
                          else
                          tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).time_spent_in_G1_per_time(time)=ceil(rand*T_G1);    
                          end
                          tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).num_of_cells_in_G1_per_time(time)=tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).num_of_cells_in_G1_per_time(time)+fraction*tumour_in_GCs(i,j,k).num_of_cells_in_G1_per_time(time);
                          
                          %the following block trasfers the cells in S phase from.... 
                          divider=(fraction*tumour_in_GCs(i,j,k).num_of_cells_in_S_per_time(time)+tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).num_of_cells_in_S_per_time(time));
                          if (divider~=0)
                          tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).time_spent_in_S_per_time(time)=round((tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).num_of_cells_in_S_per_time(time)*tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).time_spent_in_S_per_time(time)+fraction*tumour_in_GCs(i,j,k).num_of_cells_in_S_per_time(time)*tumour_in_GCs(i,j,k).time_spent_in_S_per_time(time))/divider);
                          else
                          tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).time_spent_in_S_per_time(time)=ceil(rand*T_S);    
                          end
                          tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).num_of_cells_in_S_per_time(time)=tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).num_of_cells_in_S_per_time(time)+fraction*tumour_in_GCs(i,j,k).num_of_cells_in_S_per_time(time);
                          
                          %the following block trasfers the cells in G2 phase from.... 
                          divider=(fraction*tumour_in_GCs(i,j,k).num_of_cells_in_G2_per_time(time)+tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).num_of_cells_in_G2_per_time(time));
                          if (divider~=0)
                          tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).time_spent_in_G2_per_time(time)=round((tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).num_of_cells_in_G2_per_time(time)*tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).time_spent_in_G2_per_time(time)+fraction*tumour_in_GCs(i,j,k).num_of_cells_in_G2_per_time(time)*tumour_in_GCs(i,j,k).time_spent_in_G2_per_time(time))/divider);
                          else
                          tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).time_spent_in_G2_per_time(time)=ceil(rand*T_G2);    
                          end
                          tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).num_of_cells_in_G2_per_time(time)=tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).num_of_cells_in_G2_per_time(time)+fraction*tumour_in_GCs(i,j,k).num_of_cells_in_G2_per_time(time);
                          
                          %the following block trasfers the cells in M phase from.... 
                          divider=(fraction*tumour_in_GCs(i,j,k).num_of_cells_in_M_per_time(time)+tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).num_of_cells_in_M_per_time(time));
                          if (divider~=0)
                          tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).time_spent_in_M_per_time(time)=round((tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).num_of_cells_in_M_per_time(time)*tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).time_spent_in_M_per_time(time)+fraction*tumour_in_GCs(i,j,k).num_of_cells_in_M_per_time(time)*tumour_in_GCs(i,j,k).time_spent_in_M_per_time(time))/divider);
                          else
                          tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).time_spent_in_M_per_time(time)=ceil(rand*T_M);    
                          end
                          tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).num_of_cells_in_M_per_time(time)=tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).num_of_cells_in_M_per_time(time)+fraction*tumour_in_GCs(i,j,k).num_of_cells_in_M_per_time(time);
                          
                          %the following block trasfers the cells in G0 phase from.... 
                          divider=(fraction*tumour_in_GCs(i,j,k).num_of_cells_in_G0_per_time(time)+tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).num_of_cells_in_G0_per_time(time));
                          if (divider~=0)
                          tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).time_spent_in_G0_per_time(time)=round((tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).num_of_cells_in_G0_per_time(time)*tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).time_spent_in_G0_per_time(time)+fraction*tumour_in_GCs(i,j,k).num_of_cells_in_G0_per_time(time)*tumour_in_GCs(i,j,k).time_spent_in_G0_per_time(time))/divider);
                          else
                          tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).time_spent_in_G0_per_time(time)=ceil(rand*T_G0);    
                          end
                          tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).num_of_cells_in_G0_per_time(time)=tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).num_of_cells_in_G0_per_time(time)+fraction*tumour_in_GCs(i,j,k).num_of_cells_in_G0_per_time(time);
                          
                          %the following block trasfers the cells in N phase from.... 
                          divider=(fraction*tumour_in_GCs(i,j,k).num_of_cells_in_N_per_time(time)+tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).num_of_cells_in_N_per_time(time));
                          if (divider~=0)
                          tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).time_spent_in_N_per_time(time)=round((tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).num_of_cells_in_N_per_time(time)*tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).time_spent_in_N_per_time(time)+fraction*tumour_in_GCs(i,j,k).num_of_cells_in_N_per_time(time)*tumour_in_GCs(i,j,k).time_spent_in_N_per_time(time))/divider);
                          else
                          tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).time_spent_in_N_per_time(time)=ceil(rand*T_N);    
                          end
                          tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).num_of_cells_in_N_per_time(time)=tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).num_of_cells_in_N_per_time(time)+fraction*tumour_in_GCs(i,j,k).num_of_cells_in_N_per_time(time);    
                           
                          
                          % what is done from here until the end before " end  %for neighbor=1:num_of_neighbors" : the dead cells due to radiation that are transferred to the neighbor are calculated  
                          if max(tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).dim_for_cells_radiated,tumour_in_GCs(i,j,k).dim_for_cells_radiated)==0   
                              tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).dim_for_cells_radiated=0;
                              tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).cells_radiated=[0 0];
                          else 
                           
                              case_for_neighb=zeros(max(max(tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).cells_radiated(:,2)), max(tumour_in_GCs(i,j,k).cells_radiated(:,2))),1);
                              for radiat_times_1=1:tumour_in_GCs(i,j,k).dim_for_cells_radiated
                                  if tumour_in_GCs(i,j,k).cells_radiated(radiat_times_1,2)~=0
                                   case_for_neighb(tumour_in_GCs(i,j,k).cells_radiated(radiat_times_1,2))=case_for_neighb(tumour_in_GCs(i,j,k).cells_radiated(radiat_times_1,2))+fraction*tumour_in_GCs(i,j,k).cells_radiated(radiat_times_1,1);
                                  end
                              end
                              
                              for radiat_times_2=1:tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).dim_for_cells_radiated
                                  if tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).cells_radiated(radiat_times_2,2)~=0
                                   case_for_neighb(tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).cells_radiated(radiat_times_2,2))=case_for_neighb(tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).cells_radiated(radiat_times_2,2))+tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).cells_radiated(radiat_times_2,1);
                                  end    
                              end
                              
                              case_found_for_neighb=0;
                              for radiat_case=1:max(max(tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).cells_radiated(:,2)), max(tumour_in_GCs(i,j,k).cells_radiated(:,2)))
                                  if  case_for_neighb(radiat_case,1)~=0
                                      case_found_for_neighb=case_found_for_neighb+1;
                                      tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).cells_radiated(case_found_for_neighb,1:2)=[case_for_neighb(radiat_case,1) radiat_case ];
                                  end
                              end
                              
                              [tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).dim_for_cells_radiated,ka]=size(tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).cells_radiated);
                              
                          end  %  if max(tumour_in_GCs(i+neighbors_i_s(neighbor),j+neighbors_j_s(neighbor),k+neighbors_k_s(neighbor)).dim_for_cells_radiated,tumour_in_GCs(i,j,k).dim_for_cells_radiated)==0   
                              
                       
                      end  %for neighbor=1:num_of_neighbors
                      
                      
                      direction=ceil((rand)*6);  % the direction to which the shift (due to the disappearance of the cell) will be done is randomly chosen
                      
                        
                      if direction==1 & (num_of_neighbors~=0)
                          if i<=(m-1)
                           for i_inside=i:m-1   
                                % code for sequentially shifting (from +x towards the GC) 
                                tumour_in_GCs(i_inside,j,k).num_of_cells_in_G1_per_time(time)=tumour_in_GCs(i_inside+1,j,k).num_of_cells_in_G1_per_time(time);
                                tumour_in_GCs(i_inside,j,k).num_of_cells_in_G2_per_time(time)=tumour_in_GCs(i_inside+1,j,k).num_of_cells_in_G2_per_time(time);
                                tumour_in_GCs(i_inside,j,k).num_of_cells_in_S_per_time(time)=tumour_in_GCs(i_inside+1,j,k).num_of_cells_in_S_per_time(time);
                                tumour_in_GCs(i_inside,j,k).num_of_cells_in_M_per_time(time)=tumour_in_GCs(i_inside+1,j,k).num_of_cells_in_M_per_time(time);
                                tumour_in_GCs(i_inside,j,k).num_of_cells_in_G0_per_time(time)=tumour_in_GCs(i_inside+1,j,k).num_of_cells_in_G0_per_time(time);
                                tumour_in_GCs(i_inside,j,k).num_of_cells_in_N_per_time(time)=tumour_in_GCs(i_inside+1,j,k).num_of_cells_in_N_per_time(time);
                                tumour_in_GCs(i_inside,j,k).dim_for_cells_radiated=tumour_in_GCs(i_inside+1,j,k).dim_for_cells_radiated;
                                tumour_in_GCs(i_inside,j,k).cells_radiated=tumour_in_GCs(i_inside+1,j,k).cells_radiated;
                                
                                %time_spent in phases are assigned only when the shifted cells is non-empty
                                if (tumour_in_GCs(i_inside,j,k).num_of_cells_in_G1_per_time(time)+tumour_in_GCs(i_inside,j,k).num_of_cells_in_G2_per_time(time)+tumour_in_GCs(i_inside,j,k).num_of_cells_in_S_per_time(time)+tumour_in_GCs(i_inside,j,k).num_of_cells_in_M_per_time(time)+tumour_in_GCs(i_inside,j,k).num_of_cells_in_G0_per_time(time)+tumour_in_GCs(i_inside,j,k).num_of_cells_in_N_per_time(time))~=0
                                tumour_in_GCs(i_inside,j,k).time_spent_in_G1_per_time(time)=tumour_in_GCs(i_inside+1,j,k).time_spent_in_G1_per_time(time);
                                tumour_in_GCs(i_inside,j,k).time_spent_in_G2_per_time(time)=tumour_in_GCs(i_inside+1,j,k).time_spent_in_G2_per_time(time);
                                tumour_in_GCs(i_inside,j,k).time_spent_in_S_per_time(time)=tumour_in_GCs(i_inside+1,j,k).time_spent_in_S_per_time(time);
                                tumour_in_GCs(i_inside,j,k).time_spent_in_M_per_time(time)=tumour_in_GCs(i_inside+1,j,k).time_spent_in_M_per_time(time);
                                tumour_in_GCs(i_inside,j,k).time_spent_in_G0_per_time(time)=tumour_in_GCs(i_inside+1,j,k).time_spent_in_G0_per_time(time);
                                tumour_in_GCs(i_inside,j,k).time_spent_in_N_per_time(time)=tumour_in_GCs(i_inside+1,j,k).time_spent_in_N_per_time(time);
                                end    
                                
                            end
                          end
                          tumour_in_GCs(m,j,k)=GC_out_of_cancer_final(m,j,k,m,time);    
                      end
                      
                      if direction==2 & (num_of_neighbors~=0)
                          if j<=(m-1)
                           for j_inside=j:m-1
                              
                               % code for sequentially shifting (from +y towards the GC) 
                               tumour_in_GCs(i,j_inside,k).num_of_cells_in_G1_per_time(time)=tumour_in_GCs(i,j_inside+1,k).num_of_cells_in_G1_per_time(time);
                               tumour_in_GCs(i,j_inside,k).num_of_cells_in_G2_per_time(time)=tumour_in_GCs(i,j_inside+1,k).num_of_cells_in_G2_per_time(time);
                               tumour_in_GCs(i,j_inside,k).num_of_cells_in_S_per_time(time)=tumour_in_GCs(i,j_inside+1,k).num_of_cells_in_S_per_time(time);
                               tumour_in_GCs(i,j_inside,k).num_of_cells_in_M_per_time(time)=tumour_in_GCs(i,j_inside+1,k).num_of_cells_in_M_per_time(time);
                               tumour_in_GCs(i,j_inside,k).num_of_cells_in_G0_per_time(time)=tumour_in_GCs(i,j_inside+1,k).num_of_cells_in_G0_per_time(time);
                               tumour_in_GCs(i,j_inside,k).num_of_cells_in_N_per_time(time)=tumour_in_GCs(i,j_inside+1,k).num_of_cells_in_N_per_time(time);
                               
                               tumour_in_GCs(i,j_inside,k).dim_for_cells_radiated=tumour_in_GCs(i,j_inside+1,k).dim_for_cells_radiated;
                               tumour_in_GCs(i,j_inside,k).cells_radiated=tumour_in_GCs(i,j_inside+1,k).cells_radiated;
                               
                               if (tumour_in_GCs(i,j_inside,k).num_of_cells_in_G1_per_time(time)+tumour_in_GCs(i,j_inside,k).num_of_cells_in_G2_per_time(time)+tumour_in_GCs(i,j_inside,k).num_of_cells_in_S_per_time(time)+tumour_in_GCs(i,j_inside,k).num_of_cells_in_M_per_time(time)+tumour_in_GCs(i,j_inside,k).num_of_cells_in_G0_per_time(time)+tumour_in_GCs(i,j_inside,k).num_of_cells_in_N_per_time(time))~=0
                               tumour_in_GCs(i,j_inside,k).time_spent_in_G1_per_time(time)=tumour_in_GCs(i,j_inside+1,k).time_spent_in_G1_per_time(time);
                               tumour_in_GCs(i,j_inside,k).time_spent_in_G2_per_time(time)=tumour_in_GCs(i,j_inside+1,k).time_spent_in_G2_per_time(time);
                               tumour_in_GCs(i,j_inside,k).time_spent_in_S_per_time(time)=tumour_in_GCs(i,j_inside+1,k).time_spent_in_S_per_time(time);
                               tumour_in_GCs(i,j_inside,k).time_spent_in_M_per_time(time)=tumour_in_GCs(i,j_inside+1,k).time_spent_in_M_per_time(time);
                               tumour_in_GCs(i,j_inside,k).time_spent_in_G0_per_time(time)=tumour_in_GCs(i,j_inside+1,k).time_spent_in_G0_per_time(time);
                               tumour_in_GCs(i,j_inside,k).time_spent_in_N_per_time(time)=tumour_in_GCs(i,j_inside+1,k).time_spent_in_N_per_time(time);
                               end
 
                           end
                          end
                          tumour_in_GCs(i,m,k)=GC_out_of_cancer_final(i,m,k,m,time);
                      end
                      
                      if direction==3 & (num_of_neighbors~=0)
                          if k<=(m-1)
                           for k_inside=k:m-1
                              
                               % code for sequentially shifting (from +z towards the GC) 
                               tumour_in_GCs(i,j,k_inside).num_of_cells_in_G1_per_time(time)=tumour_in_GCs(i,j,k_inside+1).num_of_cells_in_G1_per_time(time);
                               tumour_in_GCs(i,j,k_inside).num_of_cells_in_G2_per_time(time)=tumour_in_GCs(i,j,k_inside+1).num_of_cells_in_G2_per_time(time);
                               tumour_in_GCs(i,j,k_inside).num_of_cells_in_S_per_time(time)=tumour_in_GCs(i,j,k_inside+1).num_of_cells_in_S_per_time(time);
                               tumour_in_GCs(i,j,k_inside).num_of_cells_in_M_per_time(time)=tumour_in_GCs(i,j,k_inside+1).num_of_cells_in_M_per_time(time);
                               tumour_in_GCs(i,j,k_inside).num_of_cells_in_G0_per_time(time)=tumour_in_GCs(i,j,k_inside+1).num_of_cells_in_G0_per_time(time);
                               tumour_in_GCs(i,j,k_inside).num_of_cells_in_N_per_time(time)=tumour_in_GCs(i,j,k_inside+1).num_of_cells_in_N_per_time(time);
                               
                               tumour_in_GCs(i,j,k_inside).dim_for_cells_radiated=tumour_in_GCs(i,j,k_inside+1).dim_for_cells_radiated;
                               tumour_in_GCs(i,j,k_inside).cells_radiated=tumour_in_GCs(i,j,k_inside+1).cells_radiated;
                               
                               if (tumour_in_GCs(i,j,k_inside).num_of_cells_in_G1_per_time(time)+tumour_in_GCs(i,j,k_inside).num_of_cells_in_G2_per_time(time)+tumour_in_GCs(i,j,k_inside).num_of_cells_in_S_per_time(time)+tumour_in_GCs(i,j,k_inside).num_of_cells_in_M_per_time(time)+tumour_in_GCs(i,j,k_inside).num_of_cells_in_G0_per_time(time)+tumour_in_GCs(i,j,k_inside).num_of_cells_in_N_per_time(time))~=0
                               tumour_in_GCs(i,j,k_inside).time_spent_in_G1_per_time(time)=tumour_in_GCs(i,j,k_inside+1).time_spent_in_G1_per_time(time);
                               tumour_in_GCs(i,j,k_inside).time_spent_in_G2_per_time(time)=tumour_in_GCs(i,j,k_inside+1).time_spent_in_G2_per_time(time);
                               tumour_in_GCs(i,j,k_inside).time_spent_in_S_per_time(time)=tumour_in_GCs(i,j,k_inside+1).time_spent_in_S_per_time(time);
                               tumour_in_GCs(i,j,k_inside).time_spent_in_M_per_time(time)=tumour_in_GCs(i,j,k_inside+1).time_spent_in_M_per_time(time);
                               tumour_in_GCs(i,j,k_inside).time_spent_in_G0_per_time(time)=tumour_in_GCs(i,j,k_inside+1).time_spent_in_G0_per_time(time);
                               tumour_in_GCs(i,j,k_inside).time_spent_in_N_per_time(time)=tumour_in_GCs(i,j,k_inside+1).time_spent_in_N_per_time(time);
                               end
                               
                              
                           end
                        end
                          tumour_in_GCs(i,j,m)=GC_out_of_cancer_final(i,j,m,m,time);
                          
                      end
                      
                      if direction==4 & (num_of_neighbors~=0)
                          if i>=2
                            i_inside=i;
                            while (i_inside>=2)
                               % code for sequentially shifting (from -x towards the GC)  
                               tumour_in_GCs(i_inside,j,k).num_of_cells_in_G1_per_time(time)=tumour_in_GCs(i_inside-1,j,k).num_of_cells_in_G1_per_time(time);
                               tumour_in_GCs(i_inside,j,k).num_of_cells_in_G2_per_time(time)=tumour_in_GCs(i_inside-1,j,k).num_of_cells_in_G2_per_time(time);
                               tumour_in_GCs(i_inside,j,k).num_of_cells_in_S_per_time(time)=tumour_in_GCs(i_inside-1,j,k).num_of_cells_in_S_per_time(time);
                               tumour_in_GCs(i_inside,j,k).num_of_cells_in_M_per_time(time)=tumour_in_GCs(i_inside-1,j,k).num_of_cells_in_M_per_time(time);
                               tumour_in_GCs(i_inside,j,k).num_of_cells_in_G0_per_time(time)=tumour_in_GCs(i_inside-1,j,k).num_of_cells_in_G0_per_time(time);
                               tumour_in_GCs(i_inside,j,k).num_of_cells_in_N_per_time(time)=tumour_in_GCs(i_inside-1,j,k).num_of_cells_in_N_per_time(time);
                               
                               tumour_in_GCs(i_inside,j,k).dim_for_cells_radiated=tumour_in_GCs(i_inside-1,j,k).dim_for_cells_radiated;
                               tumour_in_GCs(i_inside,j,k).cells_radiated=tumour_in_GCs(i_inside-1,j,k).cells_radiated;
                               
                               if (tumour_in_GCs(i_inside,j,k).num_of_cells_in_G1_per_time(time)+tumour_in_GCs(i_inside,j,k).num_of_cells_in_G2_per_time(time)+tumour_in_GCs(i_inside,j,k).num_of_cells_in_S_per_time(time)+tumour_in_GCs(i_inside,j,k).num_of_cells_in_M_per_time(time)+tumour_in_GCs(i_inside,j,k).num_of_cells_in_G0_per_time(time)+tumour_in_GCs(i_inside,j,k).num_of_cells_in_N_per_time(time))~=0
                               tumour_in_GCs(i_inside,j,k).time_spent_in_G1_per_time(time)=tumour_in_GCs(i_inside-1,j,k).time_spent_in_G1_per_time(time);
                               tumour_in_GCs(i_inside,j,k).time_spent_in_G2_per_time(time)=tumour_in_GCs(i_inside-1,j,k).time_spent_in_G2_per_time(time);
                               tumour_in_GCs(i_inside,j,k).time_spent_in_S_per_time(time)=tumour_in_GCs(i_inside-1,j,k).time_spent_in_S_per_time(time);
                               tumour_in_GCs(i_inside,j,k).time_spent_in_M_per_time(time)=tumour_in_GCs(i_inside-1,j,k).time_spent_in_M_per_time(time);
                               tumour_in_GCs(i_inside,j,k).time_spent_in_G0_per_time(time)=tumour_in_GCs(i_inside-1,j,k).time_spent_in_G0_per_time(time);
                               tumour_in_GCs(i_inside,j,k).time_spent_in_N_per_time(time)=tumour_in_GCs(i_inside-1,j,k).time_spent_in_N_per_time(time);
                               end   
                               i_inside=i_inside-1;
                             
                            end
                         end
                          tumour_in_GCs(1,j,k)=GC_out_of_cancer_final(1,j,k,m,time);
                      end
                      
                     if direction==5 & (num_of_neighbors~=0)
                          if j>=2
                            j_inside=j;
                            while (j_inside>=2)
                               % code for sequentially shifting (from -y towards the GC)  
                               tumour_in_GCs(i,j_inside,k).num_of_cells_in_G1_per_time(time)=tumour_in_GCs(i,j_inside-1,k).num_of_cells_in_G1_per_time(time);
                               tumour_in_GCs(i,j_inside,k).num_of_cells_in_G2_per_time(time)=tumour_in_GCs(i,j_inside-1,k).num_of_cells_in_G2_per_time(time);
                               tumour_in_GCs(i,j_inside,k).num_of_cells_in_S_per_time(time)=tumour_in_GCs(i,j_inside-1,k).num_of_cells_in_S_per_time(time);
                               tumour_in_GCs(i,j_inside,k).num_of_cells_in_M_per_time(time)=tumour_in_GCs(i,j_inside-1,k).num_of_cells_in_M_per_time(time);
                               tumour_in_GCs(i,j_inside,k).num_of_cells_in_G0_per_time(time)=tumour_in_GCs(i,j_inside-1,k).num_of_cells_in_G0_per_time(time);
                               tumour_in_GCs(i,j_inside,k).num_of_cells_in_N_per_time(time)=tumour_in_GCs(i,j_inside-1,k).num_of_cells_in_N_per_time(time);
                               
                               tumour_in_GCs(i,j_inside,k).dim_for_cells_radiated=tumour_in_GCs(i,j_inside-1,k).dim_for_cells_radiated;
                               tumour_in_GCs(i,j_inside,k).cells_radiated=tumour_in_GCs(i,j_inside-1,k).cells_radiated;
                               
                               if (tumour_in_GCs(i,j_inside,k).num_of_cells_in_G1_per_time(time)+tumour_in_GCs(i,j_inside,k).num_of_cells_in_G2_per_time(time)+tumour_in_GCs(i,j_inside,k).num_of_cells_in_S_per_time(time)+tumour_in_GCs(i,j_inside,k).num_of_cells_in_M_per_time(time)+tumour_in_GCs(i,j_inside,k).num_of_cells_in_G0_per_time(time)+tumour_in_GCs(i,j_inside,k).num_of_cells_in_N_per_time(time))~=0
                               tumour_in_GCs(i,j_inside,k).time_spent_in_G1_per_time(time)=tumour_in_GCs(i,j_inside-1,k).time_spent_in_G1_per_time(time);
                               tumour_in_GCs(i,j_inside,k).time_spent_in_G2_per_time(time)=tumour_in_GCs(i,j_inside-1,k).time_spent_in_G2_per_time(time);
                               tumour_in_GCs(i,j_inside,k).time_spent_in_S_per_time(time)=tumour_in_GCs(i,j_inside-1,k).time_spent_in_S_per_time(time);
                               tumour_in_GCs(i,j_inside,k).time_spent_in_M_per_time(time)=tumour_in_GCs(i,j_inside-1,k).time_spent_in_M_per_time(time);
                               tumour_in_GCs(i,j_inside,k).time_spent_in_G0_per_time(time)=tumour_in_GCs(i,j_inside-1,k).time_spent_in_G0_per_time(time);
                               tumour_in_GCs(i,j_inside,k).time_spent_in_N_per_time(time)=tumour_in_GCs(i,j_inside-1,k).time_spent_in_N_per_time(time);
                               end   
                               j_inside=j_inside-1;
                             
                            end
                         end
                          tumour_in_GCs(i,1,k)=GC_out_of_cancer_final(i,1,k,m,time);
                      end
                      
                      if direction==6 & (num_of_neighbors~=0)
                          if k>=2
                            k_inside=k;
                            while (k_inside>=2)
                               % code for sequentially shifting (from -z towards the GC)  
                               tumour_in_GCs(i,j,k_inside).num_of_cells_in_G1_per_time(time)=tumour_in_GCs(i,j,k_inside-1).num_of_cells_in_G1_per_time(time);
                               tumour_in_GCs(i,j,k_inside).num_of_cells_in_G2_per_time(time)=tumour_in_GCs(i,j,k_inside-1).num_of_cells_in_G2_per_time(time);
                               tumour_in_GCs(i,j,k_inside).num_of_cells_in_S_per_time(time)=tumour_in_GCs(i,j,k_inside-1).num_of_cells_in_S_per_time(time);
                               tumour_in_GCs(i,j,k_inside).num_of_cells_in_M_per_time(time)=tumour_in_GCs(i,j,k_inside-1).num_of_cells_in_M_per_time(time);
                               tumour_in_GCs(i,j,k_inside).num_of_cells_in_G0_per_time(time)=tumour_in_GCs(i,j,k_inside-1).num_of_cells_in_G0_per_time(time);
                               tumour_in_GCs(i,j,k_inside).num_of_cells_in_N_per_time(time)=tumour_in_GCs(i,j,k_inside-1).num_of_cells_in_N_per_time(time);
                               
                               tumour_in_GCs(i,j,k_inside).dim_for_cells_radiated=tumour_in_GCs(i,j,k_inside-1).dim_for_cells_radiated;
                               tumour_in_GCs(i,j,k_inside).cells_radiated=tumour_in_GCs(i,j,k_inside-1).cells_radiated;
                               
                               if (tumour_in_GCs(i,j,k_inside).num_of_cells_in_G1_per_time(time)+tumour_in_GCs(i,j,k_inside).num_of_cells_in_G2_per_time(time)+tumour_in_GCs(i,j,k_inside).num_of_cells_in_S_per_time(time)+tumour_in_GCs(i,j,k_inside).num_of_cells_in_M_per_time(time)+tumour_in_GCs(i,j,k_inside).num_of_cells_in_G0_per_time(time)+tumour_in_GCs(i,j,k_inside).num_of_cells_in_N_per_time(time))~=0
                               tumour_in_GCs(i,j,k_inside).time_spent_in_G1_per_time(time)=tumour_in_GCs(i,j,k_inside-1).time_spent_in_G1_per_time(time);
                               tumour_in_GCs(i,j,k_inside).time_spent_in_G2_per_time(time)=tumour_in_GCs(i,j,k_inside-1).time_spent_in_G2_per_time(time);
                               tumour_in_GCs(i,j,k_inside).time_spent_in_S_per_time(time)=tumour_in_GCs(i,j,k_inside-1).time_spent_in_S_per_time(time);
                               tumour_in_GCs(i,j,k_inside).time_spent_in_M_per_time(time)=tumour_in_GCs(i,j,k_inside-1).time_spent_in_M_per_time(time);
                               tumour_in_GCs(i,j,k_inside).time_spent_in_G0_per_time(time)=tumour_in_GCs(i,j,k_inside-1).time_spent_in_G0_per_time(time);
                               tumour_in_GCs(i,j,k_inside).time_spent_in_N_per_time(time)=tumour_in_GCs(i,j,k_inside-1).time_spent_in_N_per_time(time);
                               end   
                               
                               k_inside=k_inside-1;
                             
                            end
                         end
                          tumour_in_GCs(i,j,1)=GC_out_of_cancer_final(i,j,1,m,time);
                      end
                      
                      
                  end  % tou if (tumour_in_GCs(i,j,k).num_of_all_cells_per_time(time)<=0.5*max_cell_density) && (tumour_in_GCs(i,j,k).num_of_all_cells_per_time(time)~=0)
                  
                      
                             
                end  % for k=1..
            end    % for j=1..   End of 3rd space scanning (for disappearing GCs that contain cells less than 0.5*max_cell_density)
        end   % for i=1..

        
        
        
         for i=1:m
           for j=1:m      % 4th space space scannning. It scans the GCs and when one containes cells more than 1.5*max_cell_density it creates a new cell in a randomly chosen direction and loads the excess cells to the new GC   
              for k=1:m
                
                  if (tumour_in_GCs(i,j,k).num_of_all_cells_per_time(time)>1.5*max_cell_density) 
                      
                      fraction=(tumour_in_GCs(i,j,k).num_of_all_cells_per_time(time)-max_cell_density)/tumour_in_GCs(i,j,k).num_of_all_cells_per_time(time); % finds the fraction of the cells that has to be loaded to the new GC
                      
                      cells_in_G1_for_new_gc= tumour_in_GCs(i,j,k).num_of_cells_in_G1_per_time(time)*fraction;
                      cells_in_S_for_new_gc= tumour_in_GCs(i,j,k).num_of_cells_in_S_per_time(time)*fraction;
                      cells_in_G2_for_new_gc= tumour_in_GCs(i,j,k).num_of_cells_in_G2_per_time(time)*fraction;
                      cells_in_M_for_new_gc= tumour_in_GCs(i,j,k).num_of_cells_in_M_per_time(time)*fraction;
                      cells_in_G0_for_new_gc= tumour_in_GCs(i,j,k).num_of_cells_in_G0_per_time(time)*fraction;
                      cells_in_N_for_new_gc= tumour_in_GCs(i,j,k).num_of_cells_in_N_per_time(time)*fraction;      % all necessary (cells in the phases, time_spent in phases, dim_for_cells_radiated, cells_radiated) for the new GC are calculated
                      time_spent_in_G1_for_new_gc=tumour_in_GCs(i,j,k).time_spent_in_G1_per_time(time);
                      time_spent_in_S_for_new_gc=tumour_in_GCs(i,j,k).time_spent_in_S_per_time(time);
                      time_spent_in_G2_for_new_gc=tumour_in_GCs(i,j,k).time_spent_in_G2_per_time(time);
                      time_spent_in_M_for_new_gc=tumour_in_GCs(i,j,k).time_spent_in_M_per_time(time);
                      time_spent_in_G0_for_new_gc=tumour_in_GCs(i,j,k).time_spent_in_G0_per_time(time);
                      time_spent_in_N_for_new_gc=tumour_in_GCs(i,j,k).time_spent_in_N_per_time(time);
                      dim_for_cells_radiated_for_new_gc=tumour_in_GCs(i,j,k).dim_for_cells_radiated;
                      cells_radiated_for_new_gc=tumour_in_GCs(i,j,k).cells_radiated;
                      cells_radiated_for_new_gc(:,1)=fraction*cells_radiated_for_new_gc(:,1);
                      
                      tumour_in_GCs(i,j,k).num_of_cells_in_G1_per_time(time)=(1-fraction)*tumour_in_GCs(i,j,k).num_of_cells_in_G1_per_time(time);
                      tumour_in_GCs(i,j,k).num_of_cells_in_S_per_time(time)=(1-fraction)*tumour_in_GCs(i,j,k).num_of_cells_in_S_per_time(time);
                      tumour_in_GCs(i,j,k).num_of_cells_in_G2_per_time(time)=(1-fraction)*tumour_in_GCs(i,j,k).num_of_cells_in_G2_per_time(time);
                      tumour_in_GCs(i,j,k).num_of_cells_in_M_per_time(time)=(1-fraction)*tumour_in_GCs(i,j,k).num_of_cells_in_M_per_time(time);    % the elements of the GC with the excess are updated so as the GC containes cells numbered equal to max_cell_density
                      tumour_in_GCs(i,j,k).num_of_cells_in_G0_per_time(time)=(1-fraction)*tumour_in_GCs(i,j,k).num_of_cells_in_G0_per_time(time);
                      tumour_in_GCs(i,j,k).num_of_cells_in_N_per_time(time)=(1-fraction)*tumour_in_GCs(i,j,k).num_of_cells_in_N_per_time(time);
                      tumour_in_GCs(i,j,k).cells_radiated(:,1)=(1-fraction)*tumour_in_GCs(i,j,k).cells_radiated(:,1);
                      
                      direction=ceil((rand)*6);
                      
                      if direction==1   
                         if i>=2
                             for i_inside=(m-1):i+2     
                                % code for shifting towards +x direction
                                tumour_in_GCs(i_inside,j,k).num_of_cells_in_G1_per_time(time)=tumour_in_GCs(i_inside-1,j,k).num_of_cells_in_G1_per_time(time);
                                tumour_in_GCs(i_inside,j,k).num_of_cells_in_G2_per_time(time)=tumour_in_GCs(i_inside-1,j,k).num_of_cells_in_G2_per_time(time);
                                tumour_in_GCs(i_inside,j,k).num_of_cells_in_S_per_time(time)=tumour_in_GCs(i_inside-1,j,k).num_of_cells_in_S_per_time(time);
                                tumour_in_GCs(i_inside,j,k).num_of_cells_in_M_per_time(time)=tumour_in_GCs(i_inside-1,j,k).num_of_cells_in_M_per_time(time);
                                tumour_in_GCs(i_inside,j,k).num_of_cells_in_G0_per_time(time)=tumour_in_GCs(i_inside-1,j,k).num_of_cells_in_G0_per_time(time);
                                tumour_in_GCs(i_inside,j,k).num_of_cells_in_N_per_time(time)=tumour_in_GCs(i_inside-1,j,k).num_of_cells_in_N_per_time(time);
                                
                                tumour_in_GCs(i_inside,j,k).dim_for_cells_radiated=tumour_in_GCs(i_inside-1,j,k).dim_for_cells_radiated;
                                tumour_in_GCs(i_inside,j,k).cells_radiated=tumour_in_GCs(i_inside-1,j,k).cells_radiated;
                                
                                if (tumour_in_GCs(i_inside,j,k).num_of_cells_in_G1_per_time(time)+tumour_in_GCs(i_inside,j,k).num_of_cells_in_G2_per_time(time)+tumour_in_GCs(i_inside,j,k).num_of_cells_in_S_per_time(time)+tumour_in_GCs(i_inside,j,k).num_of_cells_in_M_per_time(time)+tumour_in_GCs(i_inside,j,k).num_of_cells_in_G0_per_time(time)+tumour_in_GCs(i_inside,j,k).num_of_cells_in_N_per_time(time))~=0
                                tumour_in_GCs(i_inside,j,k).time_spent_in_G1_per_time(time)=tumour_in_GCs(i_inside-1,j,k).time_spent_in_G1_per_time(time);
                                tumour_in_GCs(i_inside,j,k).time_spent_in_G2_per_time(time)=tumour_in_GCs(i_inside-1,j,k).time_spent_in_G2_per_time(time);
                                tumour_in_GCs(i_inside,j,k).time_spent_in_S_per_time(time)=tumour_in_GCs(i_inside-1,j,k).time_spent_in_S_per_time(time);
                                tumour_in_GCs(i_inside,j,k).time_spent_in_M_per_time(time)=tumour_in_GCs(i_inside-1,j,k).time_spent_in_M_per_time(time);
                                tumour_in_GCs(i_inside,j,k).time_spent_in_G0_per_time(time)=tumour_in_GCs(i_inside-1,j,k).time_spent_in_G0_per_time(time);
                                tumour_in_GCs(i_inside,j,k).time_spent_in_N_per_time(time)=tumour_in_GCs(i_inside-1,j,k).time_spent_in_N_per_time(time);
                                end   
                                            
                             end
                           
                         end
                      
                      %code for creation the new GC in the gap created after shifting towards +x direction
                      tumour_in_GCs(i+1,j,k).num_of_cells_in_G1_per_time(time)=cells_in_G1_for_new_gc;
                      tumour_in_GCs(i+1,j,k).num_of_cells_in_G2_per_time(time)=cells_in_G2_for_new_gc;
                      tumour_in_GCs(i+1,j,k).num_of_cells_in_S_per_time(time)=cells_in_S_for_new_gc;
                      tumour_in_GCs(i+1,j,k).num_of_cells_in_M_per_time(time)=cells_in_M_for_new_gc;
                      tumour_in_GCs(i+1,j,k).num_of_cells_in_G0_per_time(time)=cells_in_G0_for_new_gc;
                      tumour_in_GCs(i+1,j,k).num_of_cells_in_N_per_time(time)=cells_in_N_for_new_gc;
                        
                      tumour_in_GCs(i+1,j,k).dim_for_cells_radiated=dim_for_cells_radiated_for_new_gc;
                      tumour_in_GCs(i+1,j,k).cells_radiated=cells_radiated_for_new_gc;
                      
                      tumour_in_GCs(i+1,j,k).time_spent_in_G1_per_time(time)=time_spent_in_G1_for_new_gc;
                      tumour_in_GCs(i+1,j,k).time_spent_in_G2_per_time(time)=time_spent_in_G2_for_new_gc;
                      tumour_in_GCs(i+1,j,k).time_spent_in_S_per_time(time)=time_spent_in_S_for_new_gc;
                      tumour_in_GCs(i+1,j,k).time_spent_in_M_per_time(time)=time_spent_in_M_for_new_gc;
                      tumour_in_GCs(i+1,j,k).time_spent_in_G0_per_time(time)=time_spent_in_G0_for_new_gc;
                      tumour_in_GCs(i+1,j,k).time_spent_in_N_per_time(time)=time_spent_in_N_for_new_gc;        
                                  
                     end
                      
                     if direction==2 
                         if j>=2
                             for j_inside=(m-1):j+2     
                                % code for shifting towards +y direction
                                tumour_in_GCs(i,j_inside,k).num_of_cells_in_G1_per_time(time)=tumour_in_GCs(i,j_inside-1,k).num_of_cells_in_G1_per_time(time);
                                tumour_in_GCs(i,j_inside,k).num_of_cells_in_G2_per_time(time)=tumour_in_GCs(i,j_inside-1,k).num_of_cells_in_G2_per_time(time);
                                tumour_in_GCs(i,j_inside,k).num_of_cells_in_S_per_time(time)=tumour_in_GCs(i,j_inside-1,k).num_of_cells_in_S_per_time(time);
                                tumour_in_GCs(i,j_inside,k).num_of_cells_in_M_per_time(time)=tumour_in_GCs(i,j_inside-1,k).num_of_cells_in_M_per_time(time);
                                tumour_in_GCs(i,j_inside,k).num_of_cells_in_G0_per_time(time)=tumour_in_GCs(i,j_inside-1,k).num_of_cells_in_G0_per_time(time);
                                tumour_in_GCs(i,j_inside,k).num_of_cells_in_N_per_time(time)=tumour_in_GCs(i,j_inside-1,k).num_of_cells_in_N_per_time(time);
                                
                                tumour_in_GCs(i,j_inside,k).dim_for_cells_radiated=tumour_in_GCs(i,j_inside-1,k).dim_for_cells_radiated;
                                tumour_in_GCs(i,j_inside,k).cells_radiated=tumour_in_GCs(i,j_inside-1,k).cells_radiated;
                                
                                if (tumour_in_GCs(i,j_inside,k).num_of_cells_in_G1_per_time(time)+tumour_in_GCs(i,j_inside,k).num_of_cells_in_G2_per_time(time)+tumour_in_GCs(i,j_inside,k).num_of_cells_in_S_per_time(time)+tumour_in_GCs(i,j_inside,k).num_of_cells_in_M_per_time(time)+tumour_in_GCs(i,j_inside,k).num_of_cells_in_G0_per_time(time)+tumour_in_GCs(i,j_inside,k).num_of_cells_in_N_per_time(time))~=0
                                tumour_in_GCs(i,j_inside,k).time_spent_in_G1_per_time(time)=tumour_in_GCs(i,j_inside-1,k).time_spent_in_G1_per_time(time);
                                tumour_in_GCs(i,j_inside,k).time_spent_in_G2_per_time(time)=tumour_in_GCs(i,j_inside-1,k).time_spent_in_G2_per_time(time);
                                tumour_in_GCs(i,j_inside,k).time_spent_in_S_per_time(time)=tumour_in_GCs(i,j_inside-1,k).time_spent_in_S_per_time(time);
                                tumour_in_GCs(i,j_inside,k).time_spent_in_M_per_time(time)=tumour_in_GCs(i,j_inside-1,k).time_spent_in_M_per_time(time);
                                tumour_in_GCs(i,j_inside,k).time_spent_in_G0_per_time(time)=tumour_in_GCs(i,j_inside-1,k).time_spent_in_G0_per_time(time);
                                tumour_in_GCs(i,j_inside,k).time_spent_in_N_per_time(time)=tumour_in_GCs(i,j_inside-1,k).time_spent_in_N_per_time(time);
                                end   
                                            
                             end
                           
                         end
                      
                      %code for creation the new GC in the gap created after shifting towards +y direction
                      tumour_in_GCs(i,j+1,k).num_of_cells_in_G1_per_time(time)=cells_in_G1_for_new_gc;
                      tumour_in_GCs(i,j+1,k).num_of_cells_in_G2_per_time(time)=cells_in_G2_for_new_gc;
                      tumour_in_GCs(i,j+1,k).num_of_cells_in_S_per_time(time)=cells_in_S_for_new_gc;
                      tumour_in_GCs(i,j+1,k).num_of_cells_in_M_per_time(time)=cells_in_M_for_new_gc;
                      tumour_in_GCs(i,j+1,k).num_of_cells_in_G0_per_time(time)=cells_in_G0_for_new_gc;
                      tumour_in_GCs(i,j+1,k).num_of_cells_in_N_per_time(time)=cells_in_N_for_new_gc;
                           
                      tumour_in_GCs(i,j+1,k).dim_for_cells_radiated=dim_for_cells_radiated_for_new_gc;
                      tumour_in_GCs(i,j+1,k).cells_radiated=cells_radiated_for_new_gc;
                      
                      tumour_in_GCs(i,j+1,k).time_spent_in_G1_per_time(time)=time_spent_in_G1_for_new_gc;
                      tumour_in_GCs(i,j+1,k).time_spent_in_G2_per_time(time)=time_spent_in_G2_for_new_gc;
                      tumour_in_GCs(i,j+1,k).time_spent_in_S_per_time(time)=time_spent_in_S_for_new_gc;
                      tumour_in_GCs(i,j+1,k).time_spent_in_M_per_time(time)=time_spent_in_M_for_new_gc;
                      tumour_in_GCs(i,j+1,k).time_spent_in_G0_per_time(time)=time_spent_in_G0_for_new_gc;
                      tumour_in_GCs(i,j+1,k).time_spent_in_N_per_time(time)=time_spent_in_N_for_new_gc;        
                                  
                     end
                     
                     
                     if direction==3 
                         if k>=2
                             for k_inside=(m-1):k+2     
                                 % code for shifting towards +z direction
                                tumour_in_GCs(i,j,k_inside).num_of_cells_in_G1_per_time(time)=tumour_in_GCs(i,j,k_inside-1).num_of_cells_in_G1_per_time(time);
                                tumour_in_GCs(i,j,k_inside).num_of_cells_in_G2_per_time(time)=tumour_in_GCs(i,j,k_inside-1).num_of_cells_in_G2_per_time(time);
                                tumour_in_GCs(i,j,k_inside).num_of_cells_in_S_per_time(time)=tumour_in_GCs(i,j,k_inside-1).num_of_cells_in_S_per_time(time);
                                tumour_in_GCs(i,j,k_inside).num_of_cells_in_M_per_time(time)=tumour_in_GCs(i,j,k_inside-1).num_of_cells_in_M_per_time(time);
                                tumour_in_GCs(i,j,k_inside).num_of_cells_in_G0_per_time(time)=tumour_in_GCs(i,j,k_inside-1).num_of_cells_in_G0_per_time(time);
                                tumour_in_GCs(i,j,k_inside).num_of_cells_in_N_per_time(time)=tumour_in_GCs(i,j,k_inside-1).num_of_cells_in_N_per_time(time);
                                
                                tumour_in_GCs(i,j,k_inside).dim_for_cells_radiated=tumour_in_GCs(i,j,k_inside-1).dim_for_cells_radiated;
                                tumour_in_GCs(i,j,k_inside).cells_radiated=tumour_in_GCs(i,j,k_inside-1).cells_radiated;
                                
                                if (tumour_in_GCs(i,j,k_inside).num_of_cells_in_G1_per_time(time)+tumour_in_GCs(i,j,k_inside).num_of_cells_in_G2_per_time(time)+tumour_in_GCs(i,j,k_inside).num_of_cells_in_S_per_time(time)+tumour_in_GCs(i,j,k_inside).num_of_cells_in_M_per_time(time)+tumour_in_GCs(i,j,k_inside).num_of_cells_in_G0_per_time(time)+tumour_in_GCs(i,j,k_inside).num_of_cells_in_N_per_time(time))~=0
                                tumour_in_GCs(i,j,k_inside).time_spent_in_G1_per_time(time)=tumour_in_GCs(i,j,k_inside-1).time_spent_in_G1_per_time(time);
                                tumour_in_GCs(i,j,k_inside).time_spent_in_G2_per_time(time)=tumour_in_GCs(i,j,k_inside-1).time_spent_in_G2_per_time(time);
                                tumour_in_GCs(i,j,k_inside).time_spent_in_S_per_time(time)=tumour_in_GCs(i,j,k_inside-1).time_spent_in_S_per_time(time);
                                tumour_in_GCs(i,j,k_inside).time_spent_in_M_per_time(time)=tumour_in_GCs(i,j,k_inside-1).time_spent_in_M_per_time(time);
                                tumour_in_GCs(i,j,k_inside).time_spent_in_G0_per_time(time)=tumour_in_GCs(i,j,k_inside-1).time_spent_in_G0_per_time(time);
                                tumour_in_GCs(i,j,k_inside).time_spent_in_N_per_time(time)=tumour_in_GCs(i,j,k_inside-1).time_spent_in_N_per_time(time);
                                end   
                                            
                             end
                           
                         end
                      
                      %code for creation the new GC in the gap created after shifting towards +z direction
                      tumour_in_GCs(i,j,k+1).num_of_cells_in_G1_per_time(time)=cells_in_G1_for_new_gc;
                      tumour_in_GCs(i,j,k+1).num_of_cells_in_G2_per_time(time)=cells_in_G2_for_new_gc;
                      tumour_in_GCs(i,j,k+1).num_of_cells_in_S_per_time(time)=cells_in_S_for_new_gc;
                      tumour_in_GCs(i,j,k+1).num_of_cells_in_M_per_time(time)=cells_in_M_for_new_gc;
                      tumour_in_GCs(i,j,k+1).num_of_cells_in_G0_per_time(time)=cells_in_G0_for_new_gc;
                      tumour_in_GCs(i,j,k+1).num_of_cells_in_N_per_time(time)=cells_in_N_for_new_gc;
                           
                      tumour_in_GCs(i,j,k+1).dim_for_cells_radiated=dim_for_cells_radiated_for_new_gc;
                      tumour_in_GCs(i,j,k+1).cells_radiated=cells_radiated_for_new_gc;
                      
                      
                      tumour_in_GCs(i,j,k+1).time_spent_in_G1_per_time(time)=time_spent_in_G1_for_new_gc;
                      tumour_in_GCs(i,j,k+1).time_spent_in_G2_per_time(time)=time_spent_in_G2_for_new_gc;
                      tumour_in_GCs(i,j,k+1).time_spent_in_S_per_time(time)=time_spent_in_S_for_new_gc;
                      tumour_in_GCs(i,j,k+1).time_spent_in_M_per_time(time)=time_spent_in_M_for_new_gc;
                      tumour_in_GCs(i,j,k+1).time_spent_in_G0_per_time(time)=time_spent_in_G0_for_new_gc;
                      tumour_in_GCs(i,j,k+1).time_spent_in_N_per_time(time)=time_spent_in_N_for_new_gc;        
                                  
                     end
                     
                     if direction==4 
                         if i<=m-1
                             for i_inside=2:i-2     
                                
                                % code for shifting towards -x direction
                                tumour_in_GCs(i_inside,j,k).num_of_cells_in_G1_per_time(time)=tumour_in_GCs(i_inside+1,j,k).num_of_cells_in_G1_per_time(time);
                                tumour_in_GCs(i_inside,j,k).num_of_cells_in_G2_per_time(time)=tumour_in_GCs(i_inside+1,j,k).num_of_cells_in_G2_per_time(time);
                                tumour_in_GCs(i_inside,j,k).num_of_cells_in_S_per_time(time)=tumour_in_GCs(i_inside+1,j,k).num_of_cells_in_S_per_time(time);
                                tumour_in_GCs(i_inside,j,k).num_of_cells_in_M_per_time(time)=tumour_in_GCs(i_inside+1,j,k).num_of_cells_in_M_per_time(time);
                                tumour_in_GCs(i_inside,j,k).num_of_cells_in_G0_per_time(time)=tumour_in_GCs(i_inside+1,j,k).num_of_cells_in_G0_per_time(time);
                                tumour_in_GCs(i_inside,j,k).num_of_cells_in_N_per_time(time)=tumour_in_GCs(i_inside+1,j,k).num_of_cells_in_N_per_time(time);
                                
                                tumour_in_GCs(i_inside,j,k).dim_for_cells_radiated=tumour_in_GCs(i_inside+1,j,k).dim_for_cells_radiated;
                                tumour_in_GCs(i_inside,j,k).cells_radiated=tumour_in_GCs(i_inside+1,j,k).cells_radiated;
                                
                                
                                if (tumour_in_GCs(i_inside,j,k).num_of_cells_in_G1_per_time(time)+tumour_in_GCs(i_inside,j,k).num_of_cells_in_G2_per_time(time)+tumour_in_GCs(i_inside,j,k).num_of_cells_in_S_per_time(time)+tumour_in_GCs(i_inside,j,k).num_of_cells_in_M_per_time(time)+tumour_in_GCs(i_inside,j,k).num_of_cells_in_G0_per_time(time)+tumour_in_GCs(i_inside,j,k).num_of_cells_in_N_per_time(time))~=0
                                tumour_in_GCs(i_inside,j,k).time_spent_in_G1_per_time(time)=tumour_in_GCs(i_inside+1,j,k).time_spent_in_G1_per_time(time);
                                tumour_in_GCs(i_inside,j,k).time_spent_in_G2_per_time(time)=tumour_in_GCs(i_inside+1,j,k).time_spent_in_G2_per_time(time);
                                tumour_in_GCs(i_inside,j,k).time_spent_in_S_per_time(time)=tumour_in_GCs(i_inside+1,j,k).time_spent_in_S_per_time(time);
                                tumour_in_GCs(i_inside,j,k).time_spent_in_M_per_time(time)=tumour_in_GCs(i_inside+1,j,k).time_spent_in_M_per_time(time);
                                tumour_in_GCs(i_inside,j,k).time_spent_in_G0_per_time(time)=tumour_in_GCs(i_inside+1,j,k).time_spent_in_G0_per_time(time);
                                tumour_in_GCs(i_inside,j,k).time_spent_in_N_per_time(time)=tumour_in_GCs(i_inside+1,j,k).time_spent_in_N_per_time(time);
                                end   
                                            
                             end
                           
                         end
                      
                       %code for creation of the new GC in the gap created after shifting towards -x direction   
                      tumour_in_GCs(i-1,j,k).num_of_cells_in_G1_per_time(time)=cells_in_G1_for_new_gc;
                      tumour_in_GCs(i-1,j,k).num_of_cells_in_G2_per_time(time)=cells_in_G2_for_new_gc;
                      tumour_in_GCs(i-1,j,k).num_of_cells_in_S_per_time(time)=cells_in_S_for_new_gc;
                      tumour_in_GCs(i-1,j,k).num_of_cells_in_M_per_time(time)=cells_in_M_for_new_gc;
                      tumour_in_GCs(i-1,j,k).num_of_cells_in_G0_per_time(time)=cells_in_G0_for_new_gc;
                      tumour_in_GCs(i-1,j,k).num_of_cells_in_N_per_time(time)=cells_in_N_for_new_gc;
                       
                      tumour_in_GCs(i-1,j,k).dim_for_cells_radiated=dim_for_cells_radiated_for_new_gc;
                      tumour_in_GCs(i-1,j,k).cells_radiated=cells_radiated_for_new_gc;
                      
                      
                      tumour_in_GCs(i-1,j,k).time_spent_in_G1_per_time(time)=time_spent_in_G1_for_new_gc;
                      tumour_in_GCs(i-1,j,k).time_spent_in_G2_per_time(time)=time_spent_in_G2_for_new_gc;
                      tumour_in_GCs(i-1,j,k).time_spent_in_S_per_time(time)=time_spent_in_S_for_new_gc;
                      tumour_in_GCs(i-1,j,k).time_spent_in_M_per_time(time)=time_spent_in_M_for_new_gc;
                      tumour_in_GCs(i-1,j,k).time_spent_in_G0_per_time(time)=time_spent_in_G0_for_new_gc;
                      tumour_in_GCs(i-1,j,k).time_spent_in_N_per_time(time)=time_spent_in_N_for_new_gc;        
                                  
                     end
                      
                     
                     
                     if direction==5 
                         if j<=m-1
                             for j_inside=2:j-2     
                                
                                 % code for shifting towards -y direction
                                tumour_in_GCs(i,j_inside,k).num_of_cells_in_G1_per_time(time)=tumour_in_GCs(i,j_inside+1,k).num_of_cells_in_G1_per_time(time);
                                tumour_in_GCs(i,j_inside,k).num_of_cells_in_G2_per_time(time)=tumour_in_GCs(i,j_inside+1,k).num_of_cells_in_G2_per_time(time);
                                tumour_in_GCs(i,j_inside,k).num_of_cells_in_S_per_time(time)=tumour_in_GCs(i,j_inside+1,k).num_of_cells_in_S_per_time(time);
                                tumour_in_GCs(i,j_inside,k).num_of_cells_in_M_per_time(time)=tumour_in_GCs(i,j_inside+1,k).num_of_cells_in_M_per_time(time);
                                tumour_in_GCs(i,j_inside,k).num_of_cells_in_G0_per_time(time)=tumour_in_GCs(i,j_inside+1,k).num_of_cells_in_G0_per_time(time);
                                tumour_in_GCs(i,j_inside,k).num_of_cells_in_N_per_time(time)=tumour_in_GCs(i,j_inside+1,k).num_of_cells_in_N_per_time(time);
                                
                                tumour_in_GCs(i,j_inside,k).dim_for_cells_radiated=tumour_in_GCs(i,j_inside+1,k).dim_for_cells_radiated;
                                tumour_in_GCs(i,j_inside,k).cells_radiated=tumour_in_GCs(i,j_inside+1,k).cells_radiated;
                                
                                if (tumour_in_GCs(i,j_inside,k).num_of_cells_in_G1_per_time(time)+tumour_in_GCs(i,j_inside,k).num_of_cells_in_G2_per_time(time)+tumour_in_GCs(i,j_inside,k).num_of_cells_in_S_per_time(time)+tumour_in_GCs(i,j_inside,k).num_of_cells_in_M_per_time(time)+tumour_in_GCs(i,j_inside,k).num_of_cells_in_G0_per_time(time)+tumour_in_GCs(i,j_inside,k).num_of_cells_in_N_per_time(time))~=0
                                tumour_in_GCs(i,j_inside,k).time_spent_in_G1_per_time(time)=tumour_in_GCs(i,j_inside+1,k).time_spent_in_G1_per_time(time);
                                tumour_in_GCs(i,j_inside,k).time_spent_in_G2_per_time(time)=tumour_in_GCs(i,j_inside+1,k).time_spent_in_G2_per_time(time);
                                tumour_in_GCs(i,j_inside,k).time_spent_in_S_per_time(time)=tumour_in_GCs(i,j_inside+1,k).time_spent_in_S_per_time(time);
                                tumour_in_GCs(i,j_inside,k).time_spent_in_M_per_time(time)=tumour_in_GCs(i,j_inside+1,k).time_spent_in_M_per_time(time);
                                tumour_in_GCs(i,j_inside,k).time_spent_in_G0_per_time(time)=tumour_in_GCs(i,j_inside+1,k).time_spent_in_G0_per_time(time);
                                tumour_in_GCs(i,j_inside,k).time_spent_in_N_per_time(time)=tumour_in_GCs(i,j_inside+1,k).time_spent_in_N_per_time(time);
                                end   
                                            
                             end
                           
                         end
                      
                       %code for creation the new GC in the gap created after shifting towards -y direction   
                      tumour_in_GCs(i,j-1,k).num_of_cells_in_G1_per_time(time)=cells_in_G1_for_new_gc;
                      tumour_in_GCs(i,j-1,k).num_of_cells_in_G2_per_time(time)=cells_in_G2_for_new_gc;
                      tumour_in_GCs(i,j-1,k).num_of_cells_in_S_per_time(time)=cells_in_S_for_new_gc;
                      tumour_in_GCs(i,j-1,k).num_of_cells_in_M_per_time(time)=cells_in_M_for_new_gc;
                      tumour_in_GCs(i,j-1,k).num_of_cells_in_G0_per_time(time)=cells_in_G0_for_new_gc;
                      tumour_in_GCs(i,j-1,k).num_of_cells_in_N_per_time(time)=cells_in_N_for_new_gc;
                       
                      tumour_in_GCs(i,j-1,k).dim_for_cells_radiated=dim_for_cells_radiated_for_new_gc;
                      tumour_in_GCs(i,j-1,k).cells_radiated=cells_radiated_for_new_gc;
                      
                      tumour_in_GCs(i,j-1,k).time_spent_in_G1_per_time(time)=time_spent_in_G1_for_new_gc;
                      tumour_in_GCs(i,j-1,k).time_spent_in_G2_per_time(time)=time_spent_in_G2_for_new_gc;
                      tumour_in_GCs(i,j-1,k).time_spent_in_S_per_time(time)=time_spent_in_S_for_new_gc;
                      tumour_in_GCs(i,j-1,k).time_spent_in_M_per_time(time)=time_spent_in_M_for_new_gc;
                      tumour_in_GCs(i,j-1,k).time_spent_in_G0_per_time(time)=time_spent_in_G0_for_new_gc;
                      tumour_in_GCs(i,j-1,k).time_spent_in_N_per_time(time)=time_spent_in_N_for_new_gc;        
                                  
                     end
                     
                     
                     
                     if direction==6 
                         if k<=m-1
                             for k_inside=2:k-2     
                                
                                % code for shifting towards -z direction
                                tumour_in_GCs(i,j,k_inside).num_of_cells_in_G1_per_time(time)=tumour_in_GCs(i,j,k_inside+1).num_of_cells_in_G1_per_time(time);
                                tumour_in_GCs(i,j,k_inside).num_of_cells_in_G2_per_time(time)=tumour_in_GCs(i,j,k_inside+1).num_of_cells_in_G2_per_time(time);
                                tumour_in_GCs(i,j,k_inside).num_of_cells_in_S_per_time(time)=tumour_in_GCs(i,j,k_inside+1).num_of_cells_in_S_per_time(time);
                                tumour_in_GCs(i,j,k_inside).num_of_cells_in_M_per_time(time)=tumour_in_GCs(i,j,k_inside+1).num_of_cells_in_M_per_time(time);
                                tumour_in_GCs(i,j,k_inside).num_of_cells_in_G0_per_time(time)=tumour_in_GCs(i,j,k_inside+1).num_of_cells_in_G0_per_time(time);
                                tumour_in_GCs(i,j,k_inside).num_of_cells_in_N_per_time(time)=tumour_in_GCs(i,j,k_inside+1).num_of_cells_in_N_per_time(time);
                                
                                tumour_in_GCs(i,j,k_inside).dim_for_cells_radiated=tumour_in_GCs(i,j,k_inside+1).dim_for_cells_radiated;
                                tumour_in_GCs(i,j,k_inside).cells_radiated=tumour_in_GCs(i,j,k_inside+1).cells_radiated;
                                
                                
                                if (tumour_in_GCs(i,j,k_inside).num_of_cells_in_G1_per_time(time)+tumour_in_GCs(i,j,k_inside).num_of_cells_in_G2_per_time(time)+tumour_in_GCs(i,j,k_inside).num_of_cells_in_S_per_time(time)+tumour_in_GCs(i,j,k_inside).num_of_cells_in_M_per_time(time)+tumour_in_GCs(i,j,k_inside).num_of_cells_in_G0_per_time(time)+tumour_in_GCs(i,j,k_inside).num_of_cells_in_N_per_time(time))~=0
                                tumour_in_GCs(i,j,k_inside).time_spent_in_G1_per_time(time)=tumour_in_GCs(i,j,k_inside+1).time_spent_in_G1_per_time(time);
                                tumour_in_GCs(i,j,k_inside).time_spent_in_G2_per_time(time)=tumour_in_GCs(i,j,k_inside+1).time_spent_in_G2_per_time(time);
                                tumour_in_GCs(i,j,k_inside).time_spent_in_S_per_time(time)=tumour_in_GCs(i,j,k_inside+1).time_spent_in_S_per_time(time);
                                tumour_in_GCs(i,j,k_inside).time_spent_in_M_per_time(time)=tumour_in_GCs(i,j,k_inside+1).time_spent_in_M_per_time(time);
                                tumour_in_GCs(i,j,k_inside).time_spent_in_G0_per_time(time)=tumour_in_GCs(i,j,k_inside+1).time_spent_in_G0_per_time(time);
                                tumour_in_GCs(i,j,k_inside).time_spent_in_N_per_time(time)=tumour_in_GCs(i,j,k_inside+1).time_spent_in_N_per_time(time);
                                end   
                                            
                             end
                           
                         end
                      
                      %code for creation the new GC in the gap created after shifting towards -z direction      
                      tumour_in_GCs(i,j,k-1).num_of_cells_in_G1_per_time(time)=cells_in_G1_for_new_gc;
                      tumour_in_GCs(i,j,k-1).num_of_cells_in_G2_per_time(time)=cells_in_G2_for_new_gc;
                      tumour_in_GCs(i,j,k-1).num_of_cells_in_S_per_time(time)=cells_in_S_for_new_gc;
                      tumour_in_GCs(i,j,k-1).num_of_cells_in_M_per_time(time)=cells_in_M_for_new_gc;
                      tumour_in_GCs(i,j,k-1).num_of_cells_in_G0_per_time(time)=cells_in_G0_for_new_gc;
                      tumour_in_GCs(i,j,k-1).num_of_cells_in_N_per_time(time)=cells_in_N_for_new_gc;
                       
                      tumour_in_GCs(i,j,k-1).dim_for_cells_radiated=dim_for_cells_radiated_for_new_gc;
                      tumour_in_GCs(i,j,k-1).cells_radiated=cells_radiated_for_new_gc;
                      
                      tumour_in_GCs(i,j,k-1).time_spent_in_G1_per_time(time)=time_spent_in_G1_for_new_gc;
                      tumour_in_GCs(i,j,k-1).time_spent_in_G2_per_time(time)=time_spent_in_G2_for_new_gc;
                      tumour_in_GCs(i,j,k-1).time_spent_in_S_per_time(time)=time_spent_in_S_for_new_gc;
                      tumour_in_GCs(i,j,k-1).time_spent_in_M_per_time(time)=time_spent_in_M_for_new_gc;
                      tumour_in_GCs(i,j,k-1).time_spent_in_G0_per_time(time)=time_spent_in_G0_for_new_gc;
                      tumour_in_GCs(i,j,k-1).time_spent_in_N_per_time(time)=time_spent_in_N_for_new_gc;        
                                  
                     end
                     
                     
                     
                     
                  end % for  if (tumour_in_GCs(i,j,k).num_of_all_cells_per_time(time)>1.5*max_cell_density)  
                  
              end  % for   k=1..
            end    % for j=1..     End of 4th space scanning (for creating a new GC neighboring with a GC that contain cells more than 1.5*max_cell_density and which takes the excess cells load)
        end   % for i=1..

        
        all_cells(time)=0;
        all_cells_in_G0(time)=0;
        all_cells_in_necrotic(time)=0;   %  these vectors are initialized.They will contain the number of cells (all,in G0, in necrotic, the proliferating) along time
        all_cells_in_prolif(time)=0;        
        for i=1:m
            for j=1:m       %5th time scanning for updating region for GC, the cells in necrotic, finding all cells, all_cells_in G0, etc along time
               for k=1:m
                    
                    % the element all_cells_per_time is updated
                    tumour_in_GCs(i,j,k).num_of_all_cells_per_time(time)=tumour_in_GCs(i,j,k).num_of_cells_in_G1_per_time(time)+tumour_in_GCs(i,j,k).num_of_cells_in_S_per_time(time)+tumour_in_GCs(i,j,k).num_of_cells_in_G2_per_time(time)+tumour_in_GCs(i,j,k).num_of_cells_in_M_per_time(time)+tumour_in_GCs(i,j,k).num_of_cells_in_G0_per_time(time)+tumour_in_GCs(i,j,k).num_of_cells_in_N_per_time(time);
                    for radiat_cases=1:tumour_in_GCs(i,j,k).dim_for_cells_radiated
                     tumour_in_GCs(i,j,k).num_of_all_cells_per_time(time)=tumour_in_GCs(i,j,k).num_of_all_cells_per_time(time)+tumour_in_GCs(i,j,k).cells_radiated(radiat_cases,1);
                    end
                    
                    % the elements num_of_prolif_cells_per_time, available_space_per_time are updated
                    tumour_in_GCs(i,j,k).num_of_prolif_cells_per_time(time)=tumour_in_GCs(i,j,k).num_of_cells_in_G1_per_time(time)+tumour_in_GCs(i,j,k).num_of_cells_in_S_per_time(time)+tumour_in_GCs(i,j,k).num_of_cells_in_G2_per_time(time)+tumour_in_GCs(i,j,k).num_of_cells_in_M_per_time(time);
                    tumour_in_GCs(i,j,k).available_space_per_time(time)=max(max_cell_density-tumour_in_GCs(i,j,k).num_of_all_cells_per_time(time),0);
                    
                    
                    % the elements prof_to_all, rest_to_all,dead_to_all are calculated
                    if tumour_in_GCs(i,j,k).num_of_all_cells_per_time(time)==0
                    tumour_in_GCs(i,j,k).prof_to_all(time)= 0;
                    tumour_in_GCs(i,j,k).rest_to_all(time)= 0;
                    tumour_in_GCs(i,j,k).dead_to_all(time)= 0;
                    else
                    tumour_in_GCs(i,j,k).prof_to_all(time)= tumour_in_GCs(i,j,k).num_of_prolif_cells_per_time(time)/tumour_in_GCs(i,j,k).num_of_all_cells_per_time(time);
                    tumour_in_GCs(i,j,k).rest_to_all(time)= tumour_in_GCs(i,j,k).num_of_cells_in_G0_per_time(time)/tumour_in_GCs(i,j,k).num_of_all_cells_per_time(time);
                    tumour_in_GCs(i,j,k).dead_to_all(time)= (tumour_in_GCs(i,j,k).num_of_all_cells_per_time(time)-tumour_in_GCs(i,j,k).num_of_cells_in_G1_per_time(time)-tumour_in_GCs(i,j,k).num_of_cells_in_S_per_time(time)-tumour_in_GCs(i,j,k).num_of_cells_in_G2_per_time(time)-tumour_in_GCs(i,j,k).num_of_cells_in_M_per_time(time)-tumour_in_GCs(i,j,k).num_of_cells_in_G0_per_time(time))/tumour_in_GCs(i,j,k).num_of_all_cells_per_time(time);
                    end
                    
                    % the elements num_of_cells_in_necrotic is updated
                    tumour_in_GCs(i,j,k).num_of_cells_in_necrotic_per_time(time)=tumour_in_GCs(i,j,k).dead_to_all(time)*tumour_in_GCs(i,j,k).num_of_all_cells_per_time(time);
                    
                    % the elements region_per_time is updated
                    if tumour_in_GCs(i,j,k).num_of_all_cells_per_time(time)==0
                         tumour_in_GCs(i,j,k).state_per_time(time)=0;
                    elseif tumour_in_GCs(i,j,k).dead_to_all(time)>0.9
                         tumour_in_GCs(i,j,k).state_per_time(time)=3;
                    elseif tumour_in_GCs(i,j,k).prof_to_all(time)> tumour_in_GCs(i,j,k).rest_to_all(time)
                         tumour_in_GCs(i,j,k).state_per_time(time)=1;
                    else
                         tumour_in_GCs(i,j,k).state_per_time(time)=2;
                    end
                     
                    all_cells(time)=all_cells(time)+tumour_in_GCs(i,j,k).num_of_all_cells_per_time(time);  %the cells in the current GC are added to all_cells for the current time step
                    all_cells_in_G0(time)=all_cells_in_G0(time)+tumour_in_GCs(i,j,k).num_of_cells_in_G0_per_time(time);  %as above
                    all_cells_in_necrotic(time)=all_cells_in_necrotic(time)+tumour_in_GCs(i,j,k).dead_to_all(time)*tumour_in_GCs(i,j,k).num_of_all_cells_per_time(time);  %as above                            *************THE VARIABLES HERE CAN BE USED WHEN PLOTTING all_cells, all_cells_in_G0, etc per time
                    all_cells_in_prolif(time)=all_cells_in_prolif(time)+tumour_in_GCs(i,j,k).num_of_prolif_cells_per_time(time); %as above
                     
                 end% for k=1..
            end    % for j=1..   End of 5th space scanning
         end   % for i=1...
    
        
         
         
 end %for time=1..duration..  End of time scanning
    
        
        
        
        
        
                          
                          
                          
                          
                          
                          
                          
                          
                          
                          
                          
        
  
    
 
 
 
