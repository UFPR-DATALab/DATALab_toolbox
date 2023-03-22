function [Roots, Clusters, Number] = Analysis_tree(Z,cutoff)

% Empty roots and nodes vectors
Roots = [];
Nodes = [];

% Vector length
N = max(max(Z(:,1:2)))/2+1;

% Clustering parameter
Number = N-find(Z(:,3)>=cutoff,1,'first')+1;

% Clusters
Clusters = cell(Number,1);

% Cluster number
Cluster = 0;

% Exploration of the node
ExplorationDown(Z(end,1),Z(end,3));
ExplorationDown(Z(end,2),Z(end,3));

    % Exploration of the subnodes of a node
    function ExplorationDown(Node,Dissimilarity)
        
        if ismember(Node,Nodes) || ismember(Node,Roots)
            return
        end
        
        % Adding of the current node in nodes list
        Nodes = [Nodes Node];        
        
        if Node <= N
            
            % Root
            Roots = [Roots Node];
           
            % Root whose distance is higher than limit
            [n,~]=find(Z(:,1:2)==Node);            
            if Z(n,3) >= cutoff
                Cluster = Cluster+1;
            end
            
            % Cluster
            Clusters{Cluster} = [Clusters{Cluster} Node];
            
        else
            
            % Nodes
            Node = Node-N;
            
            % Subnodes
            N1 = Z(Node,1);
            N2 = Z(Node,2);
            
            % Cluster index increment
            if Z(Node,3) < cutoff && Dissimilarity >= cutoff
                Cluster = Cluster+1;
            end
            
            % Dissimilarity of current node
            Dissimilarity = Z(Node,3);
            
            if N1 <= N && N2 <= N
                
                % Roots subnodes
                if N1 < N2
                    ExplorationDown(N1,Dissimilarity);
                    ExplorationDown(N2,Dissimilarity);
                else
                    ExplorationDown(N2,Dissimilarity);
                    ExplorationDown(N1,Dissimilarity);
                end
                
            else
                
                % Exploration of the subnodes
                ExplorationDown(N1,Dissimilarity);
                ExplorationDown(N2,Dissimilarity);
                
            end
            
        end
        
    end

end
