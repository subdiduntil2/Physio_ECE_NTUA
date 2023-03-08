function plot_tumour_lab(tumour_in_GCs,time)

% ADD your code here that plots the tumour in the 3d space, as the tumour is a time step equal to time (time can be 1,2,3,....duration)
% DO a scanning in 3d space such as   
      % for i=1:m          
       %     for j=1:m
         %        for k=1:m

%Use plot3, a command that plots in 3d space..
% For example, "plot3(xo,yo,zo,'r+')" plot a red '+' sign at the point with
% coordinates xo,yo,zo  

%The beginning has already been made for you!

[m,m,m]=size(tumour_in_GCs);

for i=1:m
    for j=1:m
        for k=1:m
                    if tumour_in_GCs(i,j,k).state_per_time(time)==0
                        hold on;
                        %print("nothing") %print nothing in the figure if the GC is out of the tumour
                    elseif tumour_in_GCs(i,j,k).state_per_time(time)==1
                        plot3(tumour_in_GCs(i,j,k).coords(1),tumour_in_GCs(i,j,k).coords(2),tumour_in_GCs(i,j,k).coords(3),'r+') % print a red plus for a GC that is in proliferating state
                        hold on; %the command hold on prevents from opening new figure with a new plot command
                    elseif  tumour_in_GCs(i,j,k).state_per_time(time)==2
                        plot3(tumour_in_GCs(i,j,k).coords(1),tumour_in_GCs(i,j,k).coords(2),tumour_in_GCs(i,j,k).coords(3),'y+') % print a yellow plus for a GC that is in proliferating state
                        hold on; %the command hold on prevents from opening new figure with a new plot command
                    elseif tumour_in_GCs(i,j,k).state_per_time(time)==3
                        plot3(tumour_in_GCs(i,j,k).coords(1),tumour_in_GCs(i,j,k).coords(2),tumour_in_GCs(i,j,k).coords(3),'m+') % print a magenta plus for a GC that is in proliferating state
                        hold on; %the command hold on prevents from opening new figure with a new plot command                                   % print a magenta plus for a GC that is in necrotic state                                    
                    end
        end
    end
end
view(3)
