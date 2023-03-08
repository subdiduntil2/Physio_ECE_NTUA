% ADD your code  to plot the number of all cells, all cells in G0, all cells in proliferation and all cells in necrotic state, from time step 4 until the end of the simulation
%The beginning has already been made for you!

plot([1:1:duration],all_cells,'b');
xlim([1, duration])
legend('Total number of cells')
xlabel('Time step')
hold on;
%TODO all_cells_in_G0 with yellow
plot([1:1:duration],all_cells_in_G0,'y');
xlim([1, duration])
legend('Total number of cells2')
xlabel('Time step')
hold on;
%TODO all_cells_in_prolif with red
plot([1:1:duration],all_cells_in_prolif,'r');
xlim([1, duration])
legend('Total number of cells3')
xlabel('Time step')
hold on;
%TODO all_cells_in_necrotic with magenta
plot([1:1:duration],all_cells_in_necrotic,'m');
xlim([1, duration])
legend('Total number of cells4')
xlabel('Time step')
hold on;
legend({'Total cells','Cells in G0','Cells in proliferation','Cells in necrosis'},'Location','northeast')
hold off;


