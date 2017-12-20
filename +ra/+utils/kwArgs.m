classdef kwArgs < handle
% 
% A very basic name/value pair parser for function inputs.
% For more formal control over function inputs (including validation), checkout Matlab's inputParser.
%
% Contact: jhadida [at] fmrib.ox.ac.uk

    properties
        CaseSensitive;
    end

    properties (SetAccess = private)
        parsed;
        accessed;
    end
    
    methods (Hidden)
        
        function name = field(self,name)
            if ~self.CaseSensitive
                name = lower(name);
            end
        end
    end
    
    methods
        
        % Constructor
        function self = kwArgs(varargin)
            self.clear();
            if nargin > 0
                self.parse(varargin{:});
            end
        end
        
        function clear(self)
            self.parsed        = struct();
            self.accessed      = {};
            self.CaseSensitive = false;
        end
        
        function copy(self,other)
            self.CaseSensitive = other.CaseSensitive;
            self.parsed        = other.parsed;
        end
        
        % report of which fields were set vs which were accessed 
        function reset_accessed(self)
            self.accessed = {};
        end
        function [ok,not_accessed] = access_report(self)
            
            all_fields   = fieldnames(self.parsed);
            not_accessed = setdiff( all_fields, self.accessed );
            ok = isempty(not_accessed);
                       
        end
        
        % convert parsed inputs as a cell of inputs that can typically be passed to a function
        function args = to_cell(self)
            args = struct2cell(self.parsed);
        end
        
        function parse(self,varargin)
            
            function x = temp(x), if iscell(x), x={x}; end; end
            
            % reset the list of accessed fields
            self.reset_accessed();
            
            % get inputs
            args = varargin;
            while iscell(args) && numel(args) == 1
                args = args{1};
            end
            
            % either a cell of key-values or a structure
            if iscell(args)
                
                if ~self.CaseSensitive
                    args(1:2:end) = cellfun( @lower, args(1:2:end), 'UniformOutput', false );
                end
                args(2:2:end) = cellfun( @temp, args(2:2:end), 'UniformOutput', false );
                self.parsed = struct(args{:});
                
            elseif isstruct(args)
                self.parse( struct2cell(args) );
                
            elseif isa(args,'ra.utils.kwArgs')
                self.copy(args);
                
            else
                error('Inputs should be either a cell of key-values or a structure.');
            end
        end
        
        % check field existence
        function yes = has(self,name)
            yes = isfield(self.parsed,self.field(name));
        end
        function yes = has_nonempty(self,name)
            name = self.field(name);
            yes  = isfield(self.parsed,name) && ~isempty(self.parsed.(name));
        end
        
        % sanitisation methods
        function sanitise(self,name,sanitise_fun)
            name = self.field(name);
            self.parsed.(name) = sanitise_fun(self.parsed.(name));
        end
        function sanitise_opt(self,name,sanitise_fun)
            if self.has(name)
                self.sanitise(name,sanitise_fun);
            end
        end
        
        % validation methods
        function validate(self,name,validate_fun)
            validate_fun( self.parsed.(self.field(name)) );
        end
        function validate_opt(self,name,validate_fun)
            if self.has(name)
                self.validate(name,validate_fun);
            end
        end
        
        % get field value or specified default
        function default = get(self,name,default)
            
            name = self.field(name);
            assert( nargin > 2 || isfield(self.parsed,name), sprintf('Required key "%s" not found.', name ));
            
            if isfield(self.parsed,name)
                default = self.parsed.(name);
                self.accessed{end+1} = name;
            else
                self.parsed.(name) = default;
            end
        end
        
        % like get, but deletes the field if it exists after returning its value
        function default = pop(self,name,default)
            
            name = self.field(name);
            assert( nargin > 2 || isfield(self.parsed,name), sprintf('Required key "%s" not found.', name ));
            
            if isfield(self.parsed,name)
                default = self.parsed.(name);
                self.parsed = rmfield(self.parsed,name);
                self.accessed{end+1} = name;
            end
        end
        
        % set field value
        function self = set(self,name,value)
            self.parsed.(self.field(name)) = value;
        end
        
        % remove field
        function self = rem(self,name)
            name = self.field(name);
            if isfield(self.parsed,name)
                self.parsed = rmfield(self.parsed,name);
            end
        end
        
    end
    
end

function c = struct2cell( s, recursive )
%
% c = struct2cell( s, recursive=false )
%
% Build a cell{ key, value } from a structure.
% Structure-arrays are not supported, but the script won't fail if the recursive flag is on, 
% and one of the values is a struct-array (it will just return the struct-array without conversion).
%
% Contact: jhadida [at] fmrib.ox.ac.uk
    
    assert( numel(s) == 1, 'Structure arrays not allowed.' );
    if nargin < 2, recursive = false; end
    
    f  = fieldnames(s);
    nf = length(f);
    c  = cell(1,2*nf);
    
    for i = 1:nf
        
        c{2*i-1} = f{i};
        c{2*i  } = s.( f{i} );
        
        if recursive && isstruct(c{2*i}) && numel(c{2*i}) == 1
            c{2*i} = struct2cell( c{2*i} );
        end
        
    end
    
end

