function plot(self)
	% Render network in 2D using Matlab 2015 graph class plotting
	G = self.graph;

	LWidths = 5*G.Edges.Weight/max(G.Edges.Weight);

	plot(G,'LineWidth',LWidths)
