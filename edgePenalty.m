function pen = edgePenalty(image,truth,edges)

    ix = find(edges == 1);
    
    XX = abs(image(ix));
    YY = abs(truth(ix));
    
    pen = corr(XX,YY,'type','Kendall');
end