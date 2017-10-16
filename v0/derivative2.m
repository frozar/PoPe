function Y = derivative2(X,dx,n,periodic)
	 
  Y=zeros(size(X));
  
  switch n
	 
    case 0
	 
      if(periodic==1)
	[n1 n2]=size(X);
	modes_n1 = [-n1/2:n1/2-1];
	W_n1 = 1i*2*(pi/(n1*dx))*modes_n1'*ones(1,n2);
	Y = real(ifft(fftshift((W_n1.^2).*fftshift(fft(X)))));
      else	  
	'use of fft in derivative2 without periodicity'
      end 
      
    case 1

      Y(2:end-1,:) = ( X(3:end,:) - 2*X(2:end-1,:) + X(1:end-2,:) ) / (dx^2);
      
      if(periodic==1)
	Y(1,:) = (X(2,:) -2*X(1,:) + X(end,:)) / (dx^2);
	Y(end,:) = (X(1,:) - 2*X(end,:) + X(end-1,:)) / (dx^2);
      else
	Y(1,:) = (X(2,:) -2*X(1,:) + X(end,:)) / (dx^2);
	Y(end,:) = (X(1,:) - 2*X(end,:) + X(end-1,:)) / (dx^2);
      end
      
    case 2
	 
      Y(3:end-2,:) = ( -X(5:end,:) + 16*X(4:end-1,:)...
		       -30*X(3:end-2,:) ...
		       +16*X(2:end-3,:) - X(1:end-4,:))/(dx^2*12);

      if(periodic==1)
	Y(1,:)     = ( -X(3,:)     + 16*X(2,:)     - 30*X(1,:)     + 16*X(end,:)     - X(end-1,:))/(dx^2*12);
	Y(2,:)     = ( -X(4,:)     + 16*X(3,:)     - 30*X(2,:)     + 16*X(1,:)       - X(end,:))/(dx^2*12);
	Y(end-1,:) = ( -X(1,:)     + 16*X(end,:)   - 30*X(end-1,:) + 16*X(end-2,:)   - X(end-3,:))/(dx^2*12);
	Y(end,:)   = ( -X(2,:)     + 16*X(1,:)     - 30*X(end,:)   + 16*X(end-1,:)   - X(end-2,:))/(dx^2*12);
      else
	Y(1,:)     = ( 11*X(5,:)   - 56*X(4,:)     +114*X(3,:)     -104*X(2,:)     +35*X(1,:)   )/(dx^2*12);
	Y(2,:)     = (   -X(5,:)   +  4*X(4,:)     +  6*X(3,:)     - 20*X(2,:)     +11*X(1,:)   )/(dx^2*12);
	Y(end-1,:) = ( 11*X(end,:) - 20*X(end-1,:) +  6*X(end-2,:) +  4*X(end-3,:) -   X(end-4,:))/(dx^2*12);
	Y(end,:)   = ( 35*X(end,:) -104*X(end-1,:) +114*X(end-2,:) - 56*X(end-3,:) +11*X(end-4,:))/(dx^2*12);
      end
      
    case 3
	 
      Y(4:end-3,:) = ( 2*X(7:end,:) - 27*X(6:end-1,:) + 270*X(5:end-2,:) ... 
		       -490*X(4:end-3,:) ...
		       +270*X(3:end-4,:) - 27*X(2:end-5,:) + 2*X(1:end-6,:))/(dx^2*180);
      
      if(periodic==1)
	Y(3,:) = ( 2*X(6,:) - 27*X(5,:) + 270*X(4,:) ... 
		   -490*X(3,:)...
		   +270*X(2,:) - 27*X(1,:) + 2*X(end,:))/(dx^2*180);
	Y(2,:) = ( 2*X(5,:) - 27*X(4,:) + 270*X(3,:) ... 
		   -490*X(2,:)...
		   +270*X(1,:) - 27*X(end,:) + 2*X(end-1,:))/(dx^2*180);
	Y(1,:) = ( 2*X(4,:) - 27*X(3,:) + 270*X(2,:) ... 
		   -490*X(1,:)...
		   +270*X(end,:) - 27*X(end-1,:) + 2*X(end-2,:))/(dx^2*180);
	Y(end,:)   = ( 2*X(3,:) - 27*X(2,:) + 270*X(1,:) ... 
		       -490*X(end,:)...
		       +270*X(end-1,:) - 27*X(end-2,:) + 2*X(end-3,:))/(dx^2*180);
	Y(end-1,:) = ( 2*X(2,:) - 27*X(1,:) + 270*X(end,:) ... 
		       -490*X(end-1,:)...
		       +270*X(end-2,:) - 27*X(end-3,:) + 2*X(end-4,:))/(dx^2*180);
	Y(end-2,:) = ( 2*X(1,:) - 27*X(end,:) + 270*X(end-1,:) ... 
		       -490*X(end-2,:)...
		       +270*X(end-3,:) - 27*X(end-4,:) + 2*X(end-5,:))/(dx^2*180);
      else      
	Y(3,:) = (    2*X(7,:)  -12*X(6,:)  +15*X(5,:) +200*X(4,:) -420*X(3,:) +228*X(2,:)  -13*X(1,:) )/(dx^2*180);
	Y(2,:) = (  -13*X(7,:)  +93*X(6,:) -285*X(5,:) +470*X(4,:) -255*X(3,:) -147*X(2,:) +137*X(1,:) )/(dx^2*180);
	Y(1,:) = (  137*X(7,:) -972*X(6,:)+2970*X(5,:)-5080*X(4,:)+5265*X(3,:)-3132*X(2,:) +812*X(1,:) )/(dx^2*180);
	Y(end-2,:) = (    2*X(end-6,:)  -12*X(end-5,:)  +15*X(end-4,:) +200*X(end-3,:) -420*X(end-2,:) +228*X(end-1,:)  -13*X(end,:) )/(dx^2*180);
	Y(end-1,:) = (  -13*X(end-6,:)  +93*X(end-5,:) -285*X(end-4,:) +470*X(end-3,:) -255*X(end-2,:) -147*X(end-1,:) +137*X(end,:) )/(dx^2*180);
	Y(end,:)   = (  137*X(end-6,:) -972*X(end-5,:)+2970*X(end-4,:)-5080*X(end-3,:)+5265*X(end-2,:)-3132*X(end-1,:) +812*X(end,:) )/(dx^2*180);
      end

  end
  
