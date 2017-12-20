classdef network < matlab.mixin.Copyable
	% Class containing information about a network of ROIs

	properties
		nodes % Table of node properties
		edges % Table of edge properties
	end

	properties(Dependent)
		n_nodes % Number of nodes/ROIs in the network
		n_edges % Number of edges
	end

	methods(Static)
		B = import(fname); % Import a saved file from the bundles folder

		function B = null(n) % n nodes, no connections
			B = ra.network(zeros(n),zeros(n));
		end

		function B = uniform(n) % n nodes, all connected to each other with distance 1 connections
			B = ra.network(ones(n)-eye(n),ones(n)-eye(n));
		end

		B = line(n_units); % Synthetic network with line connectivity

	end

	methods

		function self = network(conn,dist,shortname,longname,is_cortical,coordinates)
			% Create a network object given a connection matrix and a distance matrix
			% And optionally, cell arrays of names, and coordinates, and a logical array specifying if
			% cortical or not

			if isa(conn,'digraph')
				G = conn;
				e = [reshape(findnode(G,G.Edges.EndNodes),[],2) G.Edges.Weight G.Edges.Distance];
				self.edges = table(e(:,1),e(:,2),e(:,3),e(:,4),'VariableNames',{'From','To','Weight','Distance'});
				self.nodes = G.Nodes;
				return
			end

			n = size(conn,1);

		    if nargin < 6 || isempty(coordinates)
		        for j = 1:size(conn,1) 
					coordinates{j} = [];
		        end
		    end
		        
			if nargin < 5 || isempty(is_cortical) 
				is_cortical = logical(zeros(n,1)); % Not cortical by default
		    end

			if nargin < 3 || isempty(shortname) 
				shortname = cell(n,1);
				for j = 1:size(conn,1) 
					shortname{j} = sprintf('ROI %d',j); % Default name
				end
		    end
			
		    if nargin < 4 || isempty(longname)
		        for j = 1:size(conn,1) 
		            longname{j} = shortname{j};
		        end
		    end
		    
			self.nodes = table(shortname(:),longname(:),is_cortical,coordinates(:),'VariableNames',{'Name','Longname','Cortical','Coordinates'});

			[from,to] = meshgrid(1:n,1:n);

			self.edges = table(from(:),to(:),conn(:),dist(:),'VariableNames',{'From','To','Weight','Distance'});
			self.edges(self.edges.Weight==0,:) = [];
		end

		function [conn,dist] = netmats(self)
			conn = sparse(self.edges.From,self.edges.To,self.edges.Weight,self.n_nodes,self.n_nodes);
			dist = sparse(self.edges.From,self.edges.To,self.edges.Distance,self.n_nodes,self.n_nodes);
			conn = full(conn);
			dist = full(dist);
		end

		% strength and degree from Brain Connectivity Toolbox
		function str = strength(self)
			% Return node strengths
			CIJ = self.netmats;
			str = sum(CIJ);        % strength
		end

		function deg = degree(self)
			CIJ = self.netmats;
			CIJ = double(CIJ~=0);
			deg = sum(CIJ);
		end

		function G = digraph(self)
			% Return directed graph
			EdgeTable = table([self.edges.From,self.edges.To],self.edges.Weight,self.edges.Distance,'VariableNames',{'EndNodes','Weight','Distance'});
			G = digraph(EdgeTable,self.nodes);
		end

		function G = graph(self)
			% Return undirected graph 
			% Weights/adjacency matrix must be symmetric
			assert(unique(self.netmats-self.netmats')==0,'Network matrix is not symmetric, so cannot make an undirected graph (use net.digraph())');
			EdgeTable = table([self.edges.From,self.edges.To],self.edges.Weight,self.edges.Distance,'VariableNames',{'EndNodes','Weight','Distance'});
			EdgeTable = EdgeTable(EdgeTable.EndNodes(:,1)<=EdgeTable.EndNodes(:,2),:);
			G = graph(EdgeTable,self.nodes);
		end

		function t = named_edges(self)
			t = table(self.nodes{self.edges{:,1},1},self.nodes{self.edges{:,2},1},self.edges.Weight,self.edges.Distance,'VariableNames',{'From','To','Weight','Distance'});
		end

		function n = get.n_nodes(self)
			n = size(self.nodes,1);
		end

		function n = get.n_edges(self)
			n = size(self.edges,1);
		end

	end
end


