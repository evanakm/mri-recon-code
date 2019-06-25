function noise = createCorrelatedNoise(sigma,alpha,sizeX,sizeY,startX,startY)

    % alpha is the coefficient of autocorrelation.

    hitTop = (startX == 1);
    hitLeft = (startY == 1);
    hitBottom = (startX == sizeX);
    hitRight = (startY == sizeY);
    
    stepsToTop = startX - 1;
    stepsToBottom = sizeX - startX;
    
    stepsToLeft = startY - 1;
    stepsToRight = sizeY - startY;
    
    counter = 0;
    
    noise = zeros(sizeX,sizeY);
    noise(startX,startY) = normrnd(0,sigma);
    
%    while (~hitTop || ~hitLeft || ~hitBottom || ~hitRight)
    while (counter<=max([stepsToTop+stepsToRight,...
            stepsToTop+stepsToLeft,stepsToBottom+stepsToRight,...
            stepsToBottom+stepsToLeft]))
        
        counter = counter + 1;
        
        if(~hitTop)
            %First test to see if we've hit the top
            hitTop = (startX - counter == 1);
            currentTop = startX - counter;
        else
            currentTop = 1;
        end
        
        if(~hitBottom)
            hitBottom = (startX + counter == sizeX);
            currentBottom = startX + counter;
        else
            currentBottom = sizeX;
        end
        
        if(~hitLeft)
            hitLeft = (startX - counter == 1);
            currentTop = startX - counter;
        else
            currentTop = 1;
        end
        
        if(~hitRight)
            hitRight = (startX + counter == sizeY);
            currentRight = startY + counter;
        else
            currentRight = sizeY;
        end
        
        %Need to do this by cases, because the logic changes once we hit an
        %edge.
        %
        % In everything that follows,
        % UL = length of upper left edge
        % UR = upper right, LL = lower left, LR = lower right.
        
        if(~hitTop && ~hitLeft)
            UL = counter + 1;
            TopL = [startX - counter, startY];
        elseif(~hitTop && hitLeft)
            UL = stepsToLeft + 1;
            TopL = [startX - counter, startY];
        elseif(hitTop && ~hitLeft)
            UL = stepsToTop + 1;
            TopL = [1, startY - counter + stepsToTop];
        else
            UL = stepsToTop + stepsToLeft - counter + 1;
            TopL = [1, startY - counter + stepsToTop];
        end

        if(~hitTop && ~hitRight)
            UR = counter + 1;
            TopR = [startX - counter, startY];
        elseif(~hitTop && hitRight)
            UR = stepsToRight + 1;
            TopR = [startX - counter, startY];
        elseif(hitTop && ~hitRight)
            UR = stepsToTop + 1;
            TopR = [1, startY + counter - stepsToTop];
        else
            UR = stepsToTop + stepsToRight - counter + 1;
            TopR = [1, startY + counter - stepsToTop];
        end

        if(~hitBottom && ~hitLeft)
            LL = counter + 1;
            BottomL = [startX + counter, startY];
        elseif(~hitBottom && hitLeft)
            LL = stepsToLeft + 1;
            BottomL = [startX + counter, startY];
        elseif(hitBottom && ~hitLeft)
            LL = stepsToBottom + 1;
            BottomL = [sizeX, startY - counter + stepsToBottom];
        else
            LL = stepsToBottom + stepsToLeft - counter + 1;
            BottomL = [sizeX, startY - counter + stepsToBottom];
        end

        if(~hitBottom && ~hitRight)
            LR = counter + 1;
            BottomR = [startX + counter, startY];
        elseif(~hitBottom && hitRight)
            LR = stepsToRight + 1;
            BottomR = [startX + counter, startY];
        elseif(hitBottom && ~hitRight)
            LR = stepsToBottom + 1;
            BottomR = [sizeX, startY + counter - stepsToBottom];
        else
            LR = stepsToTop + stepsToRight - counter + 1;
            BottomR = [sizeX, startY + counter - stepsToBottom];
        end        
        
        for i=1:UL
            XX = TopL(1) + (i - 1);
            YY = TopL(2) - (i - 1);
            if(i==1 && counter <= stepsToTop)
                noise(XX,YY) = normrnd(0,sigma) + alpha*noise(XX+1,YY);
            elseif(i==UL && counter <= stepsToLeft)
                noise(XX,YY) = normrnd(0,sigma) + alpha*noise(XX,YY+1);
            else
                noise(XX,YY) = normrnd(0,sigma) + 0.5*alpha*(noise(XX+1,YY)+noise(XX,YY+1));
            end                
            
        end
        
        for i=1:UR
            XX = TopR(1) + (i - 1);
            YY = TopR(2) + (i - 1);
            if(i==1 && counter <= stepsToTop)
                noise(XX,YY) = normrnd(0,sigma) + alpha*noise(XX+1,YY);
            elseif(i==UR && counter <= stepsToRight)
                noise(XX,YY) = normrnd(0,sigma) + alpha*noise(XX,YY-1);
            else
                noise(XX,YY) = normrnd(0,sigma) + 0.5*alpha*(noise(XX+1,YY)+noise(XX,YY-1));
            end                
            
        end

        for i=1:LL
            XX = BottomL(1) - (i - 1);
            YY = BottomL(2) - (i - 1);
            if(i==1 && counter <= stepsToBottom)
                noise(XX,YY) = normrnd(0,sigma) + alpha*noise(XX-1,YY);
            elseif(i==LL && counter <= stepsToLeft)
                noise(XX,YY) = normrnd(0,sigma) + alpha*noise(XX,YY+1);
            else
                noise(XX,YY) = normrnd(0,sigma) + 0.5*alpha*(noise(XX-1,YY)+noise(XX,YY+1));
            end                
            
        end        

        for i=1:LR
            XX = BottomR(1) - (i - 1);
            YY = BottomR(2) + (i - 1);
            if(i==1 && counter <= stepsToBottom)
                noise(XX,YY) = normrnd(0,sigma) + alpha*noise(XX-1,YY);
            elseif(i==LR && counter <= stepsToLeft)
                noise(XX,YY) = normrnd(0,sigma) + alpha*noise(XX,YY-1);
            else
                noise(XX,YY) = normrnd(0,sigma) + 0.5*alpha*(noise(XX-1,YY)+noise(XX,YY-1));
            end                
            
        end        
        
        
    end
    
end