classdef (Abstract) matRad_StfGeneratorParticleRayBixelAbstract < matRad_StfGeneratorExternalRayBixelAbstract
% matRad_StfGeneratorParticleRayBixelAbstract: Abstract Superclass for
% Particle Stf Generators using the ray-bixel mechanism
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2024 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% References
%   spill structure and timing informations:
%   http://cdsweb.cern.ch/record/1182954
%   http://iopscience.iop.org/article/10.1088/0031-9155/56/20/003/meta
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    properties
        useRangeShifter = false;
        esTime = 3 * 10^6;              % [\mu s], time required for synchrotron to recharge it' spill
        spillRechargeTime = 2 * 10^6;   % [\mu s], number of particles generated in each spill
        spillSize = 4 * 10 ^ 10;        % speed of synchrotron's lateral scanning in an IES
        scanSpeed = 10;                 % [m/s], number of particles per second
        spillIntensity = 4 * 10 ^ 8;
    end

    properties (Access = protected)
        availableEnergies
        availablePeakPos
        availablePeakPosRaShi
        maxPBwidth
        pbMargin
    end

    methods
        function this = matRad_StfGeneratorParticleRayBixelAbstract(pln)
            % Constructs ExternalStfGenerator with or without pln
            if nargin < 1
                pln = [];
            end
            this@matRad_StfGeneratorExternalRayBixelAbstract(pln);

            if isempty(this.radiationMode)
                this.radiationMode = 'protons';
            end
        end

        function setDefaults(this)
            % Set default values for ExternalStfGenerator
            this.setDefaults@matRad_StfGeneratorExternalRayBixelAbstract();
        end

        function stf = writeSequencingInformation(this,stf)
            wOffset = 0;
            for f = 1:numel(this.gantryAngles)
                
                E = sort(unique([stf(f).ray(:).energy]), 'descend');
                stf(f).sequencingInfo.orderToSS = zeros(stf(f).totalNumOfBixels, 1);
                stf(f).sequencingInfo.orderToSTF = zeros(stf(f).totalNumOfBixels, 1);
                stf(f).sequencingInfo.e = zeros(stf(f).totalNumOfBixels, 1);
                % find index of each IES
                for Eix = 1:length(E)
                    s = 1;
                    stf(f).sequencingInfo.IES(Eix).energy = E(Eix);
                    for j = 1:stf(f).numOfRays
                        if(any(stf(f).ray(j).energy == E(Eix)))
                            stf(f).sequencingInfo.IES(Eix).x(s) = stf(f).ray(j).rayPos_bev(1);
                            stf(f).sequencingInfo.IES(Eix).y(s) = stf(f).ray(j).rayPos_bev(3);
                            stf(f).sequencingInfo.IES(Eix).rayIndex(s) = j;
                            stf(f).sequencingInfo.IES(Eix).wIndex(s) = wOffset +  sum(stf(f).numOfBixelsPerRay(1:(j-1))) + find(stf(f).ray(j).energy == E(Eix)); % store index
                            s = s + 1;  
                        end
                    end
                end
                % write Spot scanning order
                wOffset = wOffset + sum(stf(f).numOfBixelsPerRay);
            end
            wOffset = 0;
            for f = 1:numel(this.gantryAngles)
                E =  unique([stf(f).ray(:).energy]);
                order_count = 1;
                for Eix = 1:length(E)
                    y_sorted = sort(unique(stf(f).sequencingInfo.IES(Eix).y), 'descend');
                    x_sorted = sort(stf(f).sequencingInfo.IES(Eix).x, 'ascend');
                    for k = 1:length(y_sorted)
                        ind_y  = find(stf(f).sequencingInfo.IES(Eix).y == y_sorted(k));
                        % since backforth fasion is zig zag like, flip the order every
                        % second row
                        if ~rem(k,2)
                            ind_y   = fliplr(ind_y );
                        end
                        % according to spot scanning order, sorts w index of all
                        % bixels, use this order to transfer STF order to Spot
                        % Scanning order
                        stf(f).sequencingInfo.e(order_count:1:order_count+numel(ind_y )-1) = Eix;
                        w_ind = stf(f).sequencingInfo.IES(Eix).wIndex(ind_y) - wOffset;
                        stf(f).sequencingInfo.orderToSS(order_count:1:order_count+numel(ind_y )-1) = w_ind;
                        % according to STF order, gives us order of irradiation of
                        % each bixel, use this order to transfer Spot Scanning
                        % order to STF order
                        % orderToSTF(orderToSS) = orderToSS(orderToSTF) = 1:#bixels
                        stf(f).sequencingInfo.orderToSTF(w_ind) = order_count:1:order_count+numel(ind_y )-1;
                            
                        order_count  = order_count + numel(ind_y );
                    end
                end
                wOffset = wOffset + stf(f).totalNumOfBixels;
            end
        end

        function stf = writeSequencingTime(this,stf,w)
            steerTime = [stf(:).bixelWidth] * (10 ^ 3)/ this.scanSpeed; % [\mu s]
            spillUsage = 0;
            wOffset = 0;
            for f = 1:numel(this.gantryAngles)
                stf(f).sequencingInfo.time = zeros(stf(f).totalNumOfBixels, 1);
                E =  unique([stf(f).ray(:).energy]);
                t = 0;
                order_count = 1;
                for Eix = 1:numel(E)
                    y_sorted = sort(unique(stf(f).sequencingInfo.IES(Eix).y), 'descend');
                    x_sorted = sort(stf(f).sequencingInfo.IES(Eix).x, 'ascend');
                    for k = 1:length(y_sorted)
                        y = y_sorted(k);
                        ind_y  = find(stf(f).sequencingInfo.IES(Eix).y == y);
                        % since backforth fasion is zig zag like, flip the order every
                        % second row
                        if ~rem(k,2)
                            ind_y   = fliplr(ind_y );
                        end                        
                        for is = 1:length(ind_y)
                            s = ind_y(is);
                            x = x_sorted(s);
                            % in case there were holes inside the plan "multi"
                            % multiplies the steertime to take it into account:
                            if(k == 1 && is == 1)
                                x_prev = x;
                                y_prev = y;
                            end
                            multi = abs(x_prev - x)/stf(f).bixelWidth;
                            multi = multi + abs(y_prev - y)/stf(f).bixelWidth;
                          
                            x_prev = x;
                            y_prev = y;
                            % calculating the time:
                            numOfParticles = w(stf(f).sequencingInfo.IES(Eix).wIndex(s)-wOffset)* 10^6;
                            spillTime = numOfParticles * 10^6 /this.spillIntensity;
                            t = t + multi * steerTime(f) + spillTime;
                            % taking account of the time to recharge the spill in case the required fluence was more than spill size
                            if(spillUsage + numOfParticles > this.spillSize)
                                t = t + this.spillRechargeTime;
                                spillUsage = 0;
                            end
                            spillUsage = spillUsage  + numOfParticles;
                            stf(f).sequencingInfo.time(order_count) = t;
                            order_count  = order_count + 1;
                       end
                    end
 
                end
                 stf(f).sequencingInfo.w = w(wOffset + 1: wOffset + stf(f).totalNumOfBixels);
                 wOffset = wOffset + stf(f).totalNumOfBixels;
            end
        end

        function stf = writeStequencingPhase(this,stf,numOfPhases, motionPeriod)
            phaseTime = motionPeriod * 10 ^ 6/numOfPhases;
            for f = 1:numel(this.gantryAngles)
                realTime = phaseTime;
                stf(f).sequencingInfo.phaseMatrix = zeros(length(stf(f).sequencingInfo.time),numOfPhases);
                iPhase = 1;
                iTime = 1;
                while (iTime <= length(stf(f).sequencingInfo.time))
                    if(stf(f).sequencingInfo.time(iTime) < realTime)
                        while(iTime <= length(stf(f).sequencingInfo.time) && stf(f).sequencingInfo.time(iTime) < realTime)
                            stf(f).sequencingInfo.phaseMatrix(iTime, iPhase) = 1;
                            iTime = iTime + 1;
                        end
                    else
                        iPhase = iPhase + 1;
                        if(iPhase > numOfPhases)
                            iPhase = 1;
                        end
                        realTime = realTime + phaseTime;
                    end
                end
                stf(f).sequencingInfo.phaseMatrix = stf(f).sequencingInfo.phaseMatrix(stf(f).sequencingInfo.orderToSTF,:);
                [stf(f).sequencingInfo.phaseNum,~] = find(stf(f).sequencingInfo.phaseMatrix');
                stf(f).sequencingInfo.phaseMatrix = stf(f).sequencingInfo.phaseMatrix.*stf(f).sequencingInfo.w;
            end
        end
    end
    
    methods (Access = protected)
        function initialize(this)
            this.initialize@matRad_StfGeneratorExternalRayBixelAbstract();

            %Initialize Metadata needed for stf generators
            this.availableEnergies  = [this.machine.data.energy];
            this.availablePeakPos   = [this.machine.data.peakPos] + [this.machine.data.offset];
            availableWidths         = [this.machine.data.initFocus];
            availableWidths         = [availableWidths.SisFWHMAtIso];
            this.maxPBwidth         = max(availableWidths) / 2.355;

            matRad_cfg = MatRad_Config.instance();
            if this.useRangeShifter
                %For now only a generic range shifter is used whose thickness is
                %determined by the minimum peak width to play with
                rangeShifterEqD = round(min(this.availablePeakPos)* 1.25);
                this.availablePeakPosRaShi = this.availablePeakPos - rangeShifterEqD;

                matRad_cfg.dispWarning('Use of range shifter enabled. matRad will generate a generic range shifter with WEPL %f to enable ranges below the shortest base data entry.',rangeShifterEqD);
            end

            if sum(this.availablePeakPos<0)>0
                matRad_cfg.dispError('at least one available peak position is negative - inconsistent machine file')
            end

            %Create Water equivalent cube in ct
            this.ct = matRad_calcWaterEqD(this.ct,this.radiationMode);
        end
    end

    methods (Static)
        function [available,msg] = isAvailable(pln,machine)
            % see superclass for information

            if nargin < 2
                machine = matRad_loadMachine(pln);
            end

            % Check superclass availability
            [available,msg] = matRad_StfGeneratorExternalRayBixelAbstract.isAvailable(pln,machine);

            if ~available
                return;
            end

            available = available && isstruct(machine.data);

            available = available && all(isfield(machine.data,{'energy','peakPos','initFocus','offset'}));


            if ~available
                msg = 'Your machine file is invalid and does not contain the basic fields required for photon machines!';
            else
                msg = [];
            end
        end
    end

    methods
        function stf = generate(this, ct, cst)
            stf = generate@matRad_StfGeneratorBase(this,ct, cst);
            stf =  this.writeSequencingInformation(stf);
        end
    end

end

