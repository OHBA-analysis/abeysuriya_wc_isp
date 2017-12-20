function h = fig_letter(letter,old_letter_h,fontsize)
	if nargin < 3 || isempty(fontsize)
		fontsize = 30;
	end
	
	if nargin > 1 && ~isempty(old_letter_h)
		delete(old_letter_h);
	end
	h = text(0.03,1,letter,'Units','normalized','FontSize',fontsize,'FontWeight','bold','VerticalAlignment','top');
