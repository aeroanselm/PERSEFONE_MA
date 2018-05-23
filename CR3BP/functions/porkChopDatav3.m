function [porkChop_data] = porkChopDatav3(p1,p2,dep,arr)



mu_s=astroConstants(4);				% Sun gravitational constant
day=24*3600;						% Conversion from days to seconds

DV=zeros(length(dep),length(arr));	% DV matrix initialisation
dv1=DV; dv2=DV; TOF=DV;				% Other output matrices initialisation	

for i=1:length(dep)
	for j=1:length(arr)
		dep_date=dep(i);       		% Departure date
        	arr_date=arr(j);        % Arrival date
		
		if dep_date<arr_date
            
            if p1<11
                kep_p1=uplanet(dep_date,p1);   	% Ephemerides call for p1
            else
                kep_p1=ephNEO(dep_date,p1);
            end
            
			if p2<11
				kep_p2=uplanet(arr_date,p2);	% Ephemerides call for p2
			else
				kep_p2=ephNEO(arr_date,p2);		% Ephemerides call for p2
            end             
            
			car_p1=kep2car(kep_p1,mu_s);        % State vector for p1
    		car_p2=kep2car(kep_p2,mu_s);		% State vector for p2
        		
        	R1=car_p1(1:3);						% Position vector for p1
        	R2=car_p2(1:3);						% Position vector for p2
        		
        	tof=(arr_date-dep_date);			% Time of flight expressed in days
        	tof_sec=(arr_date-dep_date)*day;    % Time of flight expressed in seconds
        		
        	[~,~,~,~,V1,V2]=lambertMR(R1,R2,tof_sec,mu_s,0,0,0,0);	% Solving Lambert's problem
        		
        	dv1(i,j)=norm(V1-car_p1(4:6));
        	dv2(i,j)=norm(car_p2(4:6)-V2);
			TOF(i,j)=tof;
		else
        	dv1(i,j)=NaN;
        	dv2(i,j)=NaN;
			TOF(i,j)=NaN;	
		end
	end	
end

dv1=dv1';
dv2=dv2';
TOF=TOF';
DV=dv1+dv2;
C31=dv1.^2;
C32=dv2.^2;
C3=C31+C32;

porkChop_data=struct('dep',dep,'arr',arr,'TOF',TOF,'dv1',dv1,'dv2',dv2,...
                    'DV',DV,'C31',C31,'C32',C32,'C3',C3);

end
		

