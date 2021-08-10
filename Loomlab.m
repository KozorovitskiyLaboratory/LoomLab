%Written by Michael F. Priest, Northwestern University, 2021

clearvars
%Step 1: Put loom file in your matlab folder/directory (get the loom file
%by downloading it from the website
%http://www.mousebrain.org/downloads.html)

%Note:To get the actual values for the data don't run the lines directly
%below, it takes far too long data=h5read('Cortex.loom','/matrix');
%data=h5read('L6_neurons.loom','/matrix');

%Create a variable for your dataset. Yeah, I could definitely set this up
%so it just reads in the filename and sets this automatically.
%datasetName='L6_neurons.loom';
datasetName='l6_r2_cns_neurons.loom';

%Get the names of the Genes
genes=h5read(datasetName,'/row_attrs/Gene');

%To get the name of the type of cell (class, subclass)
clusterNames=h5read(datasetName,'/col_attrs/ClusterName')';
clusters=h5read(datasetName,'/col_attrs/Clusters')';

%This will give you a list of all the unique clusterID names in a single
%array called clusterID
clusterID=unique(clusterNames,'stable')';
numCellTypes=length(clusterID);
numCells=length(clusterNames);

%Enter your genes of interest that you want to count along your 3 axes here
%(I could also create a 2 axis version I guess) The patterns are called
%glut, gaba and neuromod, but could easily be called Pattern1, Pattern2,
%and Pattern3
glutPattern=["Slc17a6","Slc17a7","Slc17a8"]; %x-axis
gabaPattern=["Slc32a1","Gad1","Gad2","Slc6a5","Slc6a9"]; %y-axis
neuromodPattern=["Tph2","Th","Slc6a3","Chat","Slc5a7","Slc18a3","Pnmt","Dbh","Slc18a2","Slc6a4","Fev","Hdc","Slc18a1","Tph1","Ddc"]; %z-axis

glutIndex=find(endsWith(genes,glutPattern));
gabaIndex=find(endsWith(genes,gabaPattern));
neuromodIndex=find(endsWith(genes,neuromodPattern));

gluti=length(glutIndex);
gabai=length(gabaIndex);
neuromodi=length(neuromodIndex);

glut_data=zeros(gluti, numCells);
gaba_data=zeros(gabai, numCells);
neuromod_data=zeros(neuromodi, numCells);


for i=1:gluti
    glut_data(i,:)=h5read(datasetName,'/matrix',[1 glutIndex(i)],[numCells 1])';
end


for i=1:gabai
     gaba_data(i,:)=h5read(datasetName,'/matrix',[1 gabaIndex(i)],[numCells 1])';
end

for i=1:neuromodi
     neuromod_data(i,:)=h5read(datasetName,'/matrix',[1 neuromodIndex(i)],[numCells 1])';
end


%now initialize variables that you will put stuff into. Make the length
%match the amount of clusters in your dataset
glut_pos_cluster=zeros(numCellTypes,1);
gaba_pos_cluster=zeros(numCellTypes,1);
neuromod_pos_cluster=zeros(numCellTypes,1);
total_cluster=zeros(numCellTypes,1);

%This is not currently used but if you wanted to threshold on number of
%transcripts you could do it here (or on the data types above)
glut_sum=sum(glut_data,1);
gaba_sum=sum(gaba_data,1);
neuromod_sum=sum(neuromod_data,1);

% %This gives the pos/neg identity of each cell as gabaergic,
% glutamatergic,
%neuromodulatory
glut_pos=glut_sum>0;
gaba_pos=gaba_sum>0;
neuromod_pos=neuromod_sum>0;

%Count number of cells per cluster for each type
    for j=0:numCellTypes-1
        neuromod_pos_cluster(j+1,1)=sum(neuromod_pos(clusters(1,:)==j));
        gaba_pos_cluster(j+1,1)=sum(gaba_pos(clusters(1,:)==j));
        glut_pos_cluster(j+1,1)=sum(glut_pos(clusters(1,:)==j));
        total_cluster(j+1,1)=sum(clusters(1,:)==j);
    end

%Get the proportions of each cluster/cellType that are positive for the
%groupings you care about
proportions=zeros(numCellTypes,3);
proportions(:,1)=glut_pos_cluster./total_cluster;
proportions(:,2)=gaba_pos_cluster./total_cluster;
proportions(:,3)=neuromod_pos_cluster./total_cluster;
 
%For data presentation purposes, I wanted to exclude some neuronal
%populations. This can be done in the space below.
%find spinal cord populations
SCPops=find(startsWith(clusterID,"SC"));
%find neuroblast or neuroblast like populations
test9=find(contains(clusterID,"NBL"));
test10=find(contains(clusterID,"DETPH"));
%combine the excluded populations
excludePops=vertcat(test9,test10,SCPops);
 
 %remove excluded populations
 proportions(excludePops,:)=[];
clusterID(excludePops,:)=[];
 
%Now plot all the cell populations against your genes of interest
plot3(proportions(:,1),proportions(:,2),proportions(:,3),'ok')
 set(gcf,'color','white')
xticks([0 0.5 1])
zticks([0 0.5 1])
yticks([0 0.5 1])

%You can change labels here!
xlabel("Glutamate")
ylabel("GABA")
zlabel("ACh + Monoamine")

%I wanted to 'zoom' in on part of the graph, and label the cell populations
%of interest, so that is what is below

%This sets all axes to below 0.4. It is arbitrary and for my data set and
%you can change it to whatever you like!
idThreshold=0.4;

%Plot the 'zoomed' in figure
figure;plot3(proportions(:,1),proportions(:,2),proportions(:,3),'ok')
 set(gcf,'color','white')
 axis([0 idThreshold 0 idThreshold 0 idThreshold])
 xticks([0 idThreshold/2 idThreshold])
zticks([0 idThreshold/2 idThreshold])
yticks([0 idThreshold/2 idThreshold])
hold on
zoomPops=proportions(:,1)<idThreshold & proportions(:,2)<idThreshold & proportions(:,3)<idThreshold;
%low everything neurons
z=find(zoomPops);
for i=1:length(z)
    text(proportions(z(i),1),proportions(z(i),2),proportions(z(i),3),clusterID(z(i)),'HorizontalAlignment','right');
end

%Again, you can change these labels!
xlabel("Glutamate")
ylabel("GABA")
zlabel("ACh + Monoamine")