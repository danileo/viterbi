function y0=code(u)
    
    global deg;
    
    % code: [1 (1+D^2)/(1+D+D^2)]
    S=@(sl,ul) mod([1;0]*ul+[1 1;1 0]*sl,2);
    O=@(sl,ul) mod([1;1]*ul+[0 0;1 0]*sl,2);
    
    state=[0;0]; % initial state
    y=zeros(2,length(u)+deg); % initialize output
    
    % calculate output
    for k=1:length(u)
        y(:,k)=O(state,u(k));
        state=S(state,u(k));
    end
    
    % terminate code
    for k=length(u)+1:length(u)+deg
        %keyboard
        if isequal(state,[1;0]) || isequal(state,[0;1])
            sym=1;
        else
            sym=0;
        end
        y(:,k)=O(state,sym);
        state=S(state,sym);
    end
    
    y0=reshape(y,1,[]);
end
