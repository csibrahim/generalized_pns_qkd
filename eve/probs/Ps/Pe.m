function varargout = Pe(thetaA,thetaB,thetaE,e)

    
    varargout0 = cell(1,nargout);
    varargout1 = cell(1,nargout);

    [varargout0{:}] = Pabe(thetaA,thetaB,thetaE,0,1,e);
    [varargout1{:}] = Pabe(thetaA,thetaB,thetaE,1,1,e);
    
    varargout = cell(1,nargout);
    
    for i=1:nargout
        varargout{i} = [varargout0{i} varargout1{i}]/2;
    end

end