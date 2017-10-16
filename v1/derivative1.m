function Y = derivative1(X,dx,n,periodic)

  Y=zeros(size(X));
  
  switch n
	 
    case 0
	 
      if(periodic==1)
	[n1 n2]=size(X);
	modes_n1 = [-n1/2:n1/2-1];
	W_n1 = 1i*2*(pi/(n1*dx))*modes_n1'*ones(1,n2);
	Y = real(ifft(fftshift(W_n1.*fftshift(fft(X)))));
      else
	'use of fft in derivative1 without periodicity'	
      end 
      
    case 1
	 
      %Y(2:end,:) = (X(2:end,:)-X(1:end-1,:))/dx;
      %Y(1,:)     = (X(1,:) - X(end,:))/dx;
      
      Y(2:end-1,:) = ( X(3:end,:) - X(1:end-2,:) ) / (dx*2);
      
      if(periodic==1)
 	Y(1,:) = (X(2,:) - X(end,:)) / (dx*2);
 	Y(end,:) = (X(1,:) - X(end-1,:)) / (dx*2);
      else
 	Y(1,:)   = (-3*X(1,:)   + 4*X(2,:)     - X(3,:)) / (dx*2);
 	Y(end,:) = ( 3*X(end,:) - 4*X(end-1,:) + X(end-2,:)) / (dx*2);
      end
      
    case 2
	 
      Y(3:end-2,:) = ( -X(5:end,:) + 8*X(4:end-1,:) - 8*X(2:end-3,:) + X(1:end-4,:))/(dx*12);
      
      if(periodic==1)
	Y(1,:)     = ( -X(3,:)     + 8*X(2,:)       - 8*X(end,:)     + X(end-1,:))/(dx*12);
	Y(2,:)     = ( -X(4,:)     + 8*X(3,:)       - 8*X(1,:)       + X(end,:))/(dx*12);
	Y(end-1,:) = ( -X(1,:)     + 8*X(end,:)     - 8*X(end-2,:)   + X(end-3,:))/(dx*12);
	Y(end,:)   = ( -X(2,:)     + 8*X(1,:)       - 8*X(end-1,:)   + X(end-2,:))/(dx*12);
      else
	Y(1,:)     = ( -25*X(1,:)/12 + 4*X(2,:)  -3*X(3,:) + 4*X(4,:)/3  - X(5,:)/4 ) / dx;
	Y(2,:)     = ( -X(1,:)/4  -5*X(2,:)/6  +3*X(3,:)/2 -   X(4,:)/2  + X(5,:)/12 ) / dx;
	Y(end-1,:) = ( -X(end-4,:)/12     + X(end-3,:)/2     - 3*X(end-2,:)/2   + 5*X(end-1,:)/6 + X(end,:)/4 ) / dx;
	Y(end,:) = ( +X(end-4,:)/4    - 4*X(end-3,:)/3     +3*X(end-2,:) - 4*X(end-1,:) + 25*X(end,:)/12 ) / dx;
      end
            
    case 3
	 
      Y(4:end-3,:) = ( X(7:end,:) - 9*X(6:end-1,:) + 45*X(5:end-2,:) ... 
		       - 45*X(3:end-4,:) + 9*X(2:end-5,:) - X(1:end-6,:))/(dx*60);
      
      if(periodic==1)
	Y(3,:) = ( X(6,:) - 9*X(5,:) + 45*X(4,:) ... 
		   - 45*X(2,:) + 9*X(1,:) - X(end,:))/(dx*60);
	Y(2,:) = ( X(5,:) - 9*X(4,:) + 45*X(3,:) ... 
		   - 45*X(1,:) + 9*X(end,:) - X(end-1,:))/(dx*60);
	Y(1,:) = ( X(4,:) - 9*X(3,:) + 45*X(2,:) ... 
		   - 45*X(end,:) + 9*X(end-1,:) - X(end-2,:))/(dx*60);
	Y(end,:)   = ( X(3,:) - 9*X(2,:) + 45*X(1,:) ... 
		       - 45*X(end-1,:) + 9*X(end-2,:) - X(end-3,:))/(dx*60);
	Y(end-1,:) = ( X(2,:) - 9*X(1,:) + 45*X(end,:) ... 
		       - 45*X(end-2,:) + 9*X(end-3,:) - X(end-4,:))/(dx*60);
	Y(end-2,:) = ( X(1,:) - 9*X(end,:) + 45*X(end-1,:) ... 
		       - 45*X(end-3,:) + 9*X(end-4,:) - X(end-5,:))/(dx*60);
      else
	Y(3,:) = ( X(6,:)/30 - X(5,:)/4 + X(4,:) ... 
		   - X(3,:)/3 - X(2,:)/2 + X(1,:)/20)/dx;
	
	Y(2,:) = (-X(6,:)/20 + X(5,:)/3 - X(4,:) ... 
		  + 2*X(3,:) - 13*X(2,:)/12 - X(1,:)/5)/dx;
	
	Y(1,:) = ( X(6,:)/5 - 5*X(5,:)/4 + 10*X(4,:)/3 ... 
		   - 5*X(3,:) + 5*X(2,:) - 137*X(1,:)/60)/dx;
	
	
	Y(end-2,:) = -( X(end-5,:)/30 - X(end-4,:)/4 + X(end-3,:) ... 
			- X(end-2,:)/3 - X(end-1,:)/2 + X(end,:)/20)/dx;
	
	Y(end-1,:) = -(-X(end-5,:)/20 + X(end-4,:)/3 - X(end-3,:) ... 
		       + 2*X(end-2,:) - 13*X(end-1,:)/12 - X(end,:)/5)/dx;
	
	Y(end,:) = -( X(end-5,:)/5 - 5*X(end-4,:)/4 + 10*X(end-3,:)/3 ... 
		      - 5*X(end-2,:) + 5*X(end-1,:) - 137*X(end,:)/60)/dx;
	
      end
      
    case 4
	 
      ap4=-1/280;
      ap3=4/105;
      ap2=-1/5;
      ap1=4/5;
      %      a=
      am1=-ap1;
      am2=-ap2;
      am3=-ap3;
      am4=-ap4;
      
      Y(5:end-4,:) = ( ap4*X(9:end,:)  + ap3*X(8:end-1,:) + ap2*X(7:end-2,:) + ap1*X(6:end-3,:) ... % a*X(5:end-4,:) 
		       + am1*X(4:end-5,:)+ am2*X(3:end-6,:) + am3*X(2:end-7,:) + am4*X(1:end-8,:))/dx;
      
      if(periodic==1)
	
	ap4=-1/280;
	ap3=4/105;
	ap2=-1/5;
	ap1=4/5;
	%      a=
	am1=-ap1;
	am2=-ap2;
	am3=-ap3;
	am4=-ap4;
	
	Y(4,:) = ( ap4*X(8,:)  + ap3*X(7,:) + ap2*X(6,:) + ap1*X(5,:) ... % a*X(4,:) 
		   + am1*X(3,:)+ am2*X(2,:) + am3*X(1,:) + am4*X(end,:))/dx;
	
	Y(3,:) = ( ap4*X(7,:)  + ap3*X(6,:) + ap2*X(5,:) + ap1*X(4,:) ... % a*X(3,:) 
		   + am1*X(2,:)+ am2*X(1,:) + am3*X(end,:) + am4*X(end-1,:))/dx;
	
	Y(2,:) = ( ap4*X(6,:)  + ap3*X(5,:) + ap2*X(4,:) + ap1*X(3,:) ... % a*X(2,:) 
		   + am1*X(1,:)+ am2*X(end,:) + am3*X(end-1,:) + am4*X(end-2,:))/dx;
	
	Y(1,:) = ( ap4*X(5,:)  + ap3*X(4,:) + ap2*X(3,:) + ap1*X(2,:) ... % a*X(1,:) 
		   + am1*X(end,:)+ am2*X(end-1,:) + am3*X(end-2,:) + am4*X(end-3,:))/dx;      
	
	
	Y(end-3,:) = ( ap4*X(1,:)  + ap3*X(end,:) + ap2*X(end-1,:) + ap1*X(end-2,:) ... % a*X(end-3,:) 
		       + am1*X(end-4,:)+ am2*X(end-5,:) + am3*X(end-6,:) + am4*X(end-7,:))/dx;
	
	Y(end-2,:) = ( ap4*X(2,:)  + ap3*X(1,:) + ap2*X(end,:) + ap1*X(end-1,:) ... % a*X(end-2,:) 
		       + am1*X(end-3,:)+ am2*X(end-4,:) + am3*X(end-5,:) + am4*X(end-6,:))/dx;
	
	Y(end-1,:) = ( ap4*X(3,:)  + ap3*X(2,:) + ap2*X(1,:) + ap1*X(end,:) ... % a*X(end-1,:) 
		       + am1*X(end-2,:)+ am2*X(end-3,:) + am3*X(end-4,:) + am4*X(end-5,:))/dx;
	
	Y(end,:)   = ( ap4*X(4,:)  + ap3*X(3,:) + ap2*X(2,:) + ap1*X(1,:) ... % a*X(end,:) 
		       + am1*X(end-1,:)+ am2*X(end-2,:) + am3*X(end-3,:) + am4*X(end-4,:))/dx;
	
      else
	  
	ap4= 1/280;
	ap3=-1/28;
	ap2= 1/6;
	ap1=-1/2;
	a  = 5/4;
	am1=-9/20;
	am2=-1/2;
	am3= 1/14;
	am4=-1/168;
	
	Y(4,:) = ( ap4*X(9,:)  + ap3*X(8,:) + ap2*X(7,:) + ap1*X(6,:) + a*X(5,:) ...
		   + am1*X(4,:)+ am2*X(3,:) + am3*X(2,:) + am4*X(1,:))/dx;
	
	Y(end-3,:) = -( ap4*X(end-8,:)  + ap3*X(end-7,:) + ap2*X(end-6,:) + ap1*X(end-5,:) + a*X(end-4,:)  ...
			+ am1*X(end-3,:)+ am2*X(end-2,:) + am3*X(end-1,:) + am4*X(end,:))/dx;
	
	ap4=-1/168;
	ap3= 2/35;
	ap2=-1/4;
	ap1= 2/3;
	a  =-5/4;
	am1= 2;
	am2=-19/20;
	am3=-2/7;
	am4= 1/56;
	
	Y(3,:) = ( ap4*X(9,:)  + ap3*X(8,:) + ap2*X(7,:) + ap1*X(6,:) + a*X(5,:)  ...
		   + am1*X(4,:)+ am2*X(3,:) + am3*X(2,:) + am4*X(1,:))/dx;
	
	Y(end-2,:) = -( ap4*X(end-8,:)  + ap3*X(end-7,:) + ap2*X(end-6,:) + ap1*X(end-5,:) + a*X(end-4,:)  ...
			+ am1*X(end-3,:)+ am2*X(end-2,:) + am3*X(end-1,:) + am4*X(end,:))/dx;
	
	ap4= 1/56;
	ap3=-1/6;
	ap2= 7/10;
	ap1=-7/4;
	a  = 35/12;
	am1=-7/2;
	am2= 7/2;
	am3=-223/140;
	am4=-1/8;
	
	Y(2,:) = ( ap4*X(9,:)  + ap3*X(8,:) + ap2*X(7,:) + ap1*X(6,:) + a*X(5,:)  ...
		   + am1*X(4,:)+ am2*X(3,:) + am3*X(2,:) + am4*X(1,:))/dx;
	
	Y(end-1,:) = -( ap4*X(end-8,:)  + ap3*X(end-7,:) + ap2*X(end-6,:) + ap1*X(end-5,:) + a*X(end-4,:)  ...
			+ am1*X(end-3,:)+ am2*X(end-2,:) + am3*X(end-1,:) + am4*X(end,:))/dx;
	
	ap4=-1/8;
	ap3= 8/7;
	ap2=-14/3;
	ap1= 56/5;
	a  =-35/2;
	am1= 56/3;
	am2=-14;
	am3= 8;
	am4=-761/280;
	
	Y(1,:) = ( ap4*X(9,:)  + ap3*X(8,:) + ap2*X(7,:) + ap1*X(6,:) + a*X(5,:) ...
		   + am1*X(4,:)+ am2*X(3,:) + am3*X(2,:) + am4*X(1,:))/dx;
	
	Y(end,:) = -( ap4*X(end-8,:)  + ap3*X(end-7,:) + ap2*X(end-6,:) + ap1*X(end-5,:) + a*X(end-4,:) ...
		      + am1*X(end-3,:)+ am2*X(end-2,:) + am3*X(end-1,:) + am4*X(end,:))/dx;
	
      end
      
  end
  
