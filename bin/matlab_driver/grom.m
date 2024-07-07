%#####################################
%
%# GROM class inherits from NekROM
%# v0.0.0
%
%# Ping-Hsuan Tsai
%# 2024-07-04
%
%#####################################

classdef grom < nekrom
    % Class for GROM (Galerkin Based Reduced Order Model)
    properties
        ua
        u2a
        mu
        Re
    end
    
    methods
        % Constructor
        function obj = grom(path)
            % Call the constructor of the superclass (NekROM)
            obj@nekrom(path);
        end

        function classname = str(obj)
            classname = upper(class(obj));
        end

        function obj = get_N_dim_ops(obj, nb)
            %  Get N dimensional operators and vectors
            % : nb: number of modes
            % : returns obj: NekROM object with N dimensional operators and vectors
            obj = get_N_dim_ops@nekrom(obj, nb);
        end
        

        function obj = initialize_vars(obj, Re)
            %  Initialize other properties as needed
            % : returns obj: NeROM object with initialized variables
            obj.ua    = zeros(obj.nb+1, 1); % initialize ROM averaged velocity coefficients
            obj.u2a   = zeros(obj.nb+1, obj.nb+1); % initialize ROM averaged velocity squared coefficients
            obj.Re    = Re;
            obj.mu    = 1./obj.Re;
        end

        function [rhs] = setrhs(obj, u, method)
            % Set right hand side of the G-ROM
            % : u: ROM velocity coefficients
            % : returns rhs: right hand side of the G-ROM
            if method == "semi-implicit"
                rhs = -reshape(obj.cu*u(:,1),obj.nb,obj.nb+1)*u(:,1); % nonlinear term
                rhs = rhs-obj.mu*obj.au0(2:end,1); % viscous term of the zeroth mode
            end
        end

        function dump_rom_statistics(rom, nsteps)
            ua  = rom.ua/(nsteps);
            u2a = rom.u2a/(nsteps);

            fileID = fopen("./ua_nsteps"+nsteps+"N_"+rom.nb,'w');
            fprintf(fileID,"%24.15e\n",ua);
            fclose(fileID);
                
            fileID = fopen("./u2a_nsteps"+nsteps+"N_"+rom.nb,'w');
            fprintf(fileID,"%24.15e\n",u2a);
            fclose(fileID);
        end
    end
end