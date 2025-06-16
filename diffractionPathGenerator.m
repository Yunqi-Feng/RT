function [output, multipath] = diffractionPathGenerator(Tx, Rx, edges, CADOutput, vtx, vrx, frequency, ...
    LOS_output,pseudo_LOS_output,polarization,Jones,Ptx,antenna,enablePhase)
% diffractionPathGenerator: Finds, validates, and computes NLOS single-diffraction paths.
%
% Inputs:
%   Tx, Rx:      Transmitter and receiver positions.
%   edges:       List of all edges in the environment from generate_edges.m.
%   CADOutput:   CAD model data for obstruction checking.
%   vtx, vrx:    Velocities of Tx and Rx.
%   frequency:   Carrier frequency.
%   varargin:    Optional arguments, including qTx and qRx for antenna orientation.
%
% Outputs:
%   output:      The channel parameters for the diffracted paths (Nx22 matrix).
%   multipath:   The coordinates of the diffracted paths for visualization.

    % --- Input Parser for Optional Arguments (e.g., antenna orientation) ---
    % p = inputParser;
    % addParameter(p,'qTx',struct('center', Tx, 'angle', [0 0 0]))
    % addParameter(p,'qRx',struct('center', Rx, 'angle', [0 0 0]))
    % Add other parameters like Ptx, Jones, antenna if needed, mirroring LOSOutputGenerator
    %addParameter(p, 'Ptx', 20); % Default Ptx in dBm
    %addParameter(p, 'Jones', [1;1]); % Default Jones vector
    %addParameter(p, 'antenna', 'omni'); % Default antenna
    %addParameter(p, 'enablePhase', true); % Default enable phase
    %parse(p, varargin{:});
    % qTx = p.Results.qTx;
    % qRx = p.Results.qRx;
    % qTx = p.Parameters.qTx;
    % qRx = p.Parameters.qRx;
    %Ptx = p.Results.Ptx;
    %Jones = p.Results.Jones;
    %antenna = p.Results.antenna;
    %enablePhase = p.Results.enablePhase;
    
    output = [];
    multipath = [];
    c = 3e8;
    wavelength = c / frequency;
    if ~isempty(LOS_output)% with the LOS path
        LOS_TOA=LOS_output(8);
        LOS_P=LOS_output(9);
        LOS_AAOD=LOS_output(10);
        LOS_EAOD=LOS_output(11);
        LOS_AAOA=LOS_output(12);
        LOS_EAOA=LOS_output(13);
        LOS_PHA=LOS_output(18);
        LOS_DOP=LOS_output(20);
    end
    if ~isempty(pseudo_LOS_output)% with the LOS path
        pseudo_LOS_TOA=pseudo_LOS_output(8);
        pseudo_LOS_P=pseudo_LOS_output(9);
        pseudo_LOS_AAOD=pseudo_LOS_output(10);
        pseudo_LOS_EAOD=pseudo_LOS_output(11);
        pseudo_LOS_AAOA=pseudo_LOS_output(12);
        pseudo_LOS_EAOA=pseudo_LOS_output(13);
        pseudo_LOS_PHA=pseudo_LOS_output(18);
        pseudo_LOS_DOP=pseudo_LOS_output(20);
    end
    for i = 1:size(edges, 1)
        edge_start = edges(i, 1:3);
        edge_end = edges(i, 4:6);
        
        % --- Diffraction Point Calculation ---
        % A more accurate method than midpoint is finding the point on the edge that
        % minimizes the path length Tx-edge-Rx. This is the point of stationary phase.
        % For a line segment e1->e2, the parameter 's' for the diffraction point is:
        % s = dot(Tx-e1, e2-e1)/norm(e2-e1)^2 + dot(Rx-e1, e2-e1)/norm(e2-e1)^2, then divided by 2
        % For simplicity, we can still use the midpoint, or a more advanced method.
        % Let's stick with the midpoint for this implementation's clarity.
        diffraction_point = (edge_start + edge_end) / 2;
        
        % --- Path Validation ---
        % Check for obstruction from Tx to diffraction point
        is_path_to_edge_clear = verifyPath(Tx, diffraction_point, diffraction_point - Tx, [0,0,0,0], [0,0,0,0], CADOutput, 2, false);
        
        % Check for obstruction from diffraction point to Rx
        is_path_from_edge_clear = verifyPath(diffraction_point, Rx, Rx - diffraction_point, [0,0,0,0], [0,0,0,0], CADOutput, 2, false);
        
        if is_path_to_edge_clear && is_path_from_edge_clear
            
            % --- Path Geometry and Loss Calculation ---
            d1 = norm(diffraction_point - Tx);
            d2 = norm(Rx - diffraction_point);
            distance = d1 + d2;
            delayLOS = distance / c;
            
            % Calculate diffraction loss using KED (assuming grazing incidence, h=0)
            diffraction_loss_dB = calculate_ked_loss(0, d1, d2, wavelength);
            
            % --- Start Filling the 22-column Output Vector ---
            output1 = nan(1, 22);
            output1(1) = -1; % Use -1 to identify this as a purely diffracted path
            
            % --- DoD, DoA, and Angle Calculation (AoD/AoA) ---
            dodNoRot = diffraction_point - Tx;
            doaNoRot = Rx - diffraction_point;
            
            % Apply antenna frame rotation
            % dod = coordinateRotation(dodNoRot, [0 0 0], qTx.angle, 'frame');
            % doa = coordinateRotation(doaNoRot, [0 0 0], qRx.angle, 'frame');
            dod = dodNoRot;
            doa = doaNoRot;
            
            output1(2:4) = dod; % Direction of Departure (vector)
            output1(5:7) = doa; % Direction of Arrival (vector)
            output1(8) = delayLOS;
            % [aodAz, aodEl] = vector2angle(dod);
            % [aoaAz, aoaEl] = vector2angle(doa);
            % 
            % output1(10) = aodAz; % AoD Azimuth (degrees)
            % output1(11) = aodEl; % AoD Elevation (degrees)
            % output1(12) = aoaAz; % AoA Azimuth (degrees)
            % output1(13) = aoaEl; % AoA Elevation (degrees)


            output1(10)=atan2(dod(2),dod(1))*180/pi;
            % Aod elevation
            output1(11)=atan2(dod(3),sqrt(dod(1)^2+dod(2)^2))*180/pi;
            % Aoa azimuth
            output1(12)=atan2(doa(2),doa(1))*180/pi;
            % Aoa elevation
            output1(13)=atan2(doa(3),sqrt(doa(1)^2+doa(2)^2))*180/pi;

            

            if ~isempty(LOS_output)% still assume LOS path exists
                 rel_aaod=abs(output1(10)-LOS_AAOD);
                 rel_eaod=abs(output1(11)-LOS_EAOD);
                 rel_aaoa=abs(output1(12)-LOS_AAOA);
                 rel_eaoa=abs(output1(13)-LOS_EAOA);
                 if strcmp(antenna, 'dir')
                     Gt=horn_gain(rel_eaod,rel_aaod);
                     Gr=horn_gain(rel_eaoa,rel_aaoa);
                 elseif strcmp(antenna, 'omni')
                     Gt=3.3*(cos(0.5*pi*cosd(rel_eaod))/sind(rel_eaod))^2;
                     Gr=3.3*(cos(0.5*pi*cosd(rel_eaoa))/sind(rel_eaoa))^2;
                 elseif strcmp(antenna, 'scan')
                     Gt=23.7;
                     Gr=23.7;
                 end
             elseif isempty(LOS_output) && (strcmp(antenna, 'dir') || strcmp(antenna, 'scan')) 
                 Gt=23.7;
                 Gr=23.7;
 % %If LOS path does not exist in directional channel 
 % sounding schemes, the Tx/Rx antennas will point on the 
 % best NLOS direction, but it's done in post processing
                 % LOS_AAOD=output1(10);
                 % LOS_EAOD=output1(11);
                 % LOS_AAOA=output1(12);
                 % LOS_EAOA=output1(13);
                 % LOS_output=[LOS_AAOD,LOS_EAOD,LOS_AAOA,LOS_EAOA];
             elseif isempty(LOS_output) && strcmp(antenna, 'omni')
                 rel_aaod=abs(output1(10)-pseudo_LOS_AAOD);
                 rel_eaod=abs(output1(11)-pseudo_LOS_EAOD);
                 rel_aaoa=abs(output1(12)-pseudo_LOS_AAOA);
                 rel_eaoa=abs(output1(13)-pseudo_LOS_EAOA);
                 Gt=3.3*(cos(0.5*pi*cosd(rel_eaod))/sind(rel_eaod))^2;
                 Gr=3.3*(cos(0.5*pi*cosd(rel_eaoa))/sind(rel_eaoa))^2;
             end





            
            % --- Doppler Shift Calculation ---
            % Doppler is caused by change in path length.
            % Change rate = proj. of v_tx on Tx->D + proj. of v_rx on D->Rx
            d_path_dt = dot(vtx, (Tx - diffraction_point)/d1) + dot(vrx, (diffraction_point - Rx)/d2);
            dopplerFactor = -d_path_dt / c;
            output1(20) = dopplerFactor * frequency; % Doppler Shift (Hz)

            % --- Phase Calculation ---
            if enablePhase
                output1(18) = mod(distance / wavelength * 2 * pi, 2 * pi); % Phase Shift
            else
                output1(18) = 0;
            end
            
            % --- Path Gain and Polarization ---
            % In this simplified NLOS case, assume omni antenna gain for the path finding
            %Gt = 1; Gr = 1; % Linear gain
            
            free_space_path_loss_dB = 20 * log10(4 * pi * distance / wavelength);
            total_path_loss_dB = free_space_path_loss_dB + diffraction_loss_dB;
            


            if strcmp(polarization, 'V-V') || strcmp(polarization, 'H-H')
                output1(9) = Ptx + Gt + Gr - total_path_loss_dB;
                output1(16)=output1(9)-Ptx-Gt-Gr;
                Jones_dB = 10*log10(Jones);
                output1(17) = NaN;%cross-polarized channel gain,reflection
                output1(19) = min(Jones_dB + output1(9));
                output1(22) = NaN;
            %elseif strcmp(polarization, 'cross') || strcmp(polarization, 'dual')
            elseif strcmp(polarization, 'dual')
                %XPD = 34;
                Jones_dB = 10*log10(Jones);
                output1(9)=Ptx+Gt+Gr- total_path_loss_dB;%co-polarized channel gain
                output1(16)=output1(9)-Ptx-Gt-Gr;
                output1(17)=Ptx+Gt+Gr- total_path_loss_dB;%another co-polarization channel gain
                output1(19)= min(Jones_dB(:,1) + output1(9));% cross-polarization channel gain at Tx
                output1(22)= min(Jones_dB(:,1) + output1(17)); % cross-polarization channel gain at Rx
            end

            if output1(19)<=-120
                output1(19)=NaN;
            end
            if output1(22)<=-120
                output1(22)=NaN;
            end
            if enablePhase
                output1(18) = mod(distance/wavelength*2*pi,2*pi);
            else
                output1(18) = 0;
            end
            output1(14)=Gt;
            output1(15)=Gr;
            output1(20) = dopplerFactor * frequency;


            if output1(9)<=-120 && (output1(17)<=-120 || isnan(output1(17)))
                output1=[];
                outputQd(indexOutput).output1=output1;
                multipath(iterateNumberOfRowsArraysOfPlanes,:)=[];
            else
                if output1(17)<=-120
                    output1(17)=NaN;
                end
                if output1(9)<=-120
                    output1(9)=NaN;
                end
                %outputQd(indexOutput).output1 = output1;
            end



            %output1(9) = Ptx + Gt + Gr - total_path_loss_dB; % Path Gain (dBm)
            
            % Fill other columns for consistency
            % output1(8) = delayLOS; % Time Delay (s)
            % % output1(14) = 10*log10(Gt); % Gt (dBi)
            % % output1(15) = 10*log10(Gr); % Gr (dBi)
            % output1(16) = -total_path_loss_dB; % Path Loss (dB)
            % output1(17) = NaN; % Cross-pol path gain 1
            % output1(19) = NaN; % Cross-pol path gain 2
            % output1(21) = -1;  % Cluster Index (optional, -1 for non-QD)
            % output1(22) = NaN; % Other XPol gain

            % Append the valid path to the output
            output = [output; output1];
            multipath = [multipath; [Tx, diffraction_point, Rx]];
        end
    end
end