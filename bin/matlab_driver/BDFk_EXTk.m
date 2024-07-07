    classdef BDFk_EXTk
        properties
            u
            alphas
            betas
            dt
            T_final
            nsteps
            iostep
            ito
            ext
            hfac
            rhs
        end
        methods 
            % Constructor
            function obj = BDFk_EXTk(nsteps, dt, iostep, nb)
                obj.alphas  = zeros(3,3);     % initialize vectors to hold EXTk coefficients
                obj.betas   = zeros(4,3);     % initialize vectors to hold BDFk coefficients
                obj.ext     = zeros(nb,3);    % initialize vectors to hold extrapolation vectors
                obj.u       = zeros(nb+1, 3); % initialize vectors to hold ROM coefficients
                obj.hfac    = [];
                obj.dt      = dt; % Simultaion time step size
                obj.nsteps  = nsteps; % Simulation total time steps
                obj.iostep  = iostep;  % Every #steps to store ROM quantities
                obj.T_final = obj.nsteps*obj.dt;  % Simulation final time
            end
            function obj = setup(obj)
                % Setup BDFk/EXTk coefficients

                % Setup EXT1 coefficients
                obj.alphas(1,1) =  1.0;

                % Setup EXT2 coefficients
                obj.alphas(1,2) =  2.0;
                obj.alphas(2,2) = -1.0;

                % Setup EXT3 coefficients
                obj.alphas(1,3) =  3.0;
                obj.alphas(2,3) = -3.0;
                obj.alphas(3,3) =  1.0;

                % BDF1 coefficients 
                obj.betas(1,1)  =  1.0;
                obj.betas(2,1)  = -1.0;

                % BDF2 coefficients 
                obj.betas(1,2)  =  1.5;
                obj.betas(2,2)  = -2.0;
                obj.betas(3,2)  =  0.5;

                % BDF3 coefficients 
                obj.betas(1,3)  =  11.0/6;
                obj.betas(2,3)  = -3.0;
                obj.betas(3,3)  =  1.5;
                obj.betas(4,3)  = -1.0/3;
            end

            function obj = setrhs(obj, rom, u)

                obj.ext(:,3) = obj.ext(:,2);
                obj.ext(:,2) = obj.ext(:,1);
                obj.ext(:,1) = rom.setrhs(u, "semi-implicit");

                obj.rhs = obj.ext*obj.alphas(:,obj.ito);
                % Later move rom.bu into rom.setrhs (need inverse though)
                obj.rhs = obj.rhs - rom.bu*(u(2:end,:)*obj.betas(2:end,obj.ito))/obj.dt;
            end
            
            function [next] = advance(obj, rom)
                if isempty(obj.hfac)
                    h=rom.bu*obj.betas(1,obj.ito)/obj.dt+rom.mu*rom.au;
                    obj.hfac=chol(h);
                end
                next = [1,(obj.hfac\(obj.hfac'\obj.rhs))'];
                if obj.ito <= 2
                    obj.hfac = [];
                end
            end

        end

        methods (Static)
            function a = shift(a,b,n)
                for i=n:-1:2
                    a(:,i)=a(:,i-1);
                end
                a(:,1)=b;
            end

            function rom = collect_statistics(rom, u_new)
                rom.ua  = rom.ua + u_new';
                rom.u2a = rom.u2a + u_new'*u_new;
            end
        end
    end