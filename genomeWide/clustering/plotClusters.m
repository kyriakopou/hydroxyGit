function plotClusters(T_cluster, T_out, dataDesc)

numOfClusters = length(unique(T_cluster.Cluster));

for k=1:numOfClusters
        rowsCluster = (T_cluster.Cluster == k);
        
        subplot(2, ceil(numOfClusters / 2), k);
        
        if strcmp(dataDesc, 'eff') || strcmp(dataDesc, 'both')
            %define colors
            rgbColor1 = [0.8500   0.3250    0.0980;
                        0         0.4470    0.7410;
                        0.9290    0.6940    0.1250];
                    
%             rgbColor1 = [0         0.4470    0.7410];
            
            
            effAllCluster = {[T_out(rowsCluster,:).maint_d0, T_out(rowsCluster,:).maint_d3, T_out(rowsCluster,:).maint_d6]; 
            [T_out(rowsCluster,:).deNovo_d0, T_out(rowsCluster,:).deNovo_d3, T_out(rowsCluster,:).deNovo_d6];
            [T_out(rowsCluster,:).hydroxy_d0, T_out(rowsCluster,:).hydroxy_d3, T_out(rowsCluster,:).hydroxy_d6]};
            

            %ONLY HYDROXY
%             effAllCluster = {[T_final(rowsCluster,:).hydroxy_d0, T_final(rowsCluster,:).hydroxy_d3, T_final(rowsCluster,:).hydroxy_d6]};          
%             aboxplot(effAllCluster, 'Colormap', rgbColor1, 'labels', {'d0', 'd3', 'd6'});
            
            %ONLY DE-NOVO
%             effAllCluster = {[T_final(rowsCluster,:).deNovo_d0, T_final(rowsCluster,:).deNovo_d3, T_final(rowsCluster,:).deNovo_d6]};

            aboxplot(effAllCluster, 'Colormap', rgbColor1, 'labels', {'d0', 'd3', 'd6'});
            
        elseif strcmp(dataDesc, 'levels')
            %define colors
            rgbColor2 = [0.8500    0.3250    0.0980;
                        0.929      0.694     0.125;
                        0.466      0.674     0.1880;
                        0          0.447     0.7410];
                        
            levelsAllCluster = {[T_out(rowsCluster,:).mm_d0, T_out(rowsCluster,:).mm_d3, T_out(rowsCluster,:).mm_d6]; 
            [T_out(rowsCluster,:).toth_d0, T_out(rowsCluster,:).toth_d3, T_out(rowsCluster,:).toth_d6];
            [T_out(rowsCluster,:).um_d0, T_out(rowsCluster,:).um_d3, T_out(rowsCluster,:).um_d6]
            [T_out(rowsCluster,:).uu_d0, T_out(rowsCluster,:).uu_d3, T_out(rowsCluster,:).uu_d6]};
        
            aboxplot(levelsAllCluster, 'Colormap', rgbColor2, 'labels', {'d0', 'd3', 'd6'});
        %if the clustering is made on both eff and levels
        %plot first the eff and in another figure the levels
               
        end        
              
        ylim([0,1]);
              
        %display numOfCpGs in each cluster
        numOfCpGsInCluster = sum(rowsCluster, 1);
%         avgSamplesPerCpG = sum(numOfObsPerCpG(rowsCluster)) / numOfCpGsInCluster;
        fprintf('There are %d CpGs in cluster %d \n', numOfCpGsInCluster, k);
        
end


if strcmp(dataDesc, 'both')
figure();
    for k=1:numOfClusters
        
        rowsCluster = T_cluster.Cluster == k;

        %define colors
        rgbColor2 = [0.8500   0.3250    0.0980;
                    0         0.4470    0.7410;
                    0.466     0.674     0.1880;
                    0         0.447     0.7410];


        levelsAllCluster = {[T_out(rowsCluster,:).mm_d0, T_out(rowsCluster,:).mm_d3, T_out(rowsCluster,:).mm_d6]; 
        [T_out(rowsCluster,:).toth_d0, T_out(rowsCluster,:).toth_d3, T_out(rowsCluster,:).toth_d6];
        [T_out(rowsCluster,:).um_d0, T_out(rowsCluster,:).um_d3, T_out(rowsCluster,:).um_d6]
        [T_out(rowsCluster,:).uu_d0, T_out(rowsCluster,:).uu_d3, T_out(rowsCluster,:).uu_d6]};

        aboxplot(levelsAllCluster, 'Colormap', rgbColor2, 'labels', {'d0', 'd3', 'd6'});
        
        
    end    
end


end