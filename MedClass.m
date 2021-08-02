% =========================================================================
% SUSCEPTIBILITY CLASS
% returns susceptibility of specific material in atmospheric pressure
% Output methods are:
% sDriver - returns susceptibility of neutral atoms for IR
% sXuv - returns susceptibility of neutral atoms for XUV
% sPlasma - returns susceptibility of plasma
% aLength - returns absorption length of gas for specific wavelength

classdef MedClass
    properties
        sellmeier
        kingston
        planck = 6.6260693e-34; %[J*s]
        epsilonZero = 8.854187817e-12; %[F/m]
        electronCharge = 1.60217662e-19; %[C] 
        electronMass = 9.10938356e-31; %[kg]
        speedLight = 299792458; %[m/s]
        numberDensity = 2.504e25; %[m^-3] amount of particles of ideal gas in m^3, 1 atm and 20°C
        gasConstant = 8.3144598 %[J*mol^-1*K^-1]
        temperature = 293.15; %[K], 20°C
        avogadroConstant = 6.022140857e23; %[mol^-1]
        p_a = 101325; %Atmospheric pressure [pa]
        k=1.380649*10^(-23); %Boltzmann constant [J/K]
        T=293.15; %Termodynamic temperature [K]
    end
    methods
        function obj = MedClass
            obj.sellmeier = readtable('Medium/MatData.txt','Delimiter','tab');
            obj.kingston = readtable('Medium/MatDataKingston.txt','Delimiter','tab');
        end
        function susc = sDriver(obj, gas, lambda)
            % SUSCEPTIBILITY OF NEUTRALS FOR DRIVER
            % sellmeier equation from refractiveindex.info
            % inputs
            % gas - string - "Ar", "He", "Ne", "Kr", "Xe"
            % lambda - wavelength [m]
            mask = ismember(obj.sellmeier{:,1}, gas);
            formula = str2func(strjoin(obj.sellmeier{mask,2}));
            index = formula(lambda.*1e6);   %lambda needs to be in micrometers
            susc = ((index).^2) - 1;
        end
        function susc = sDriverKing(obj, gas, lambda)
            % SUSCEPTIBILITY OF NEUTRALS FOR DRIVER
            % equations from Kingston - 1960 - the refractive indeces and constants of inert gases
            % inputs
            % gas - string - "Ar", "He", "Ne", "Kr", "Xe"
            % lambda - wavelength [m]
            mask = ismember(obj.kingston{:,1}, gas);
            formula = str2func(strjoin(obj.kingston{mask,2}));
            susc = formula(lambda.*1e10);   %lambda needs to be in angs
        end
        function susc = sXuv(obj, gas, lambda)
            % SUSCEPTIBILITY OF NEUTRALS FOR XUV
            % Calculated from formula which uses scattering factor
            % (Attwood formula 3.9)
            % Scattering factor saved in txt files
            % inputs
            % gas - string - "Ar", "He", "Ne", "Kr", "Xe"
            % lambda - wavelength [m]
            scatter = dlmread(['Medium/f1_'+gas+'.txt']);
            lambdaEnergy = (obj.planck.*obj.speedLight)./(lambda.*obj.electronCharge); % in eV
            fZero = interp1(scatter(:,1),scatter(:,2),(lambdaEnergy.*1e-3)); % energy in keV
            classicalRadius = (obj.electronCharge.^2)./(4.*pi.*obj.epsilonZero.*obj.electronMass.*(obj.speedLight.^2));
            susc = -((obj.numberDensity.*classicalRadius.*(lambda.^2))./(pi)).*fZero;
        end
        function susc = sPlasma(obj, eta, lambda)
            % SUSCEPTIBILITY OF PLASMA
            % value calculated as -(wp/w0)
            % inputs
            % gas - string - "Ar", "He", "Ne", "Kr", "Xe"
            % lambda - wavelength [m]
            pFreq = (eta.*obj.numberDensity.*obj.electronCharge.^2)./(obj.electronMass.*obj.epsilonZero);
            lFreq = ((2.*pi.*obj.speedLight)./lambda);
            susc = -(pFreq./(lFreq.^2));
        end
        function length = aLength(obj, gas, lambda, pressure)
            % inputs
            % gas - string - "Ar", "He", "Ne", "Kr", "Xe"
            % lambda - wavelength [m]
            % pressure in bar
            lambdaEnergy = (obj.planck*obj.speedLight)/(lambda*obj.electronCharge); % in eV
            
            % use scatterning factor f2 from cxro, for lower photon energy from NIST
            % as they have a little bit more data in that direction
            if lambdaEnergy > 15.9
                scatter = dlmread(['Medium/f2_'+gas+'.txt']); 
                fTwo = interp1(scatter(:,1),scatter(:,3),(lambdaEnergy));
            else
                scatter = dlmread(['Medium/f2_'+gas+'_NIST.txt']); 
                fTwo = interp1(scatter(:,1),scatter(:,2),(lambdaEnergy/1000));
            end
            
            classicalRadius = (obj.electronCharge.^2)./(4.*pi.*obj.epsilonZero.*obj.electronMass.*(obj.speedLight.^2));
            sigma = 2*classicalRadius*(lambda)*fTwo;
            particleNumber = (pressure.*1e5)./(obj.gasConstant*obj.temperature); %[mol*m^-3], pressure needs to be in [Pa]
            particleNumber = obj.avogadroConstant.*particleNumber; %[m^-3]
            length = 1/(particleNumber*sigma);
        end
    end
end