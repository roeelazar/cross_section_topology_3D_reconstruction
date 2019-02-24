function out = xor2(x)

if isempty(x)
    out = 0 ;
else
    
    out = sign( sum( xor( x(:,1) , x(:,2) ) ) ) ;
    
end

end