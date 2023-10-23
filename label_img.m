function label_img

g = uicontrol('Style','popupmenu','String', possible_labels,'Position',[20 20 100 20]);

g.Callback = @selection;




    function selection(src, event)
        val = g.Value;
        str = g.String;
        str{val};
        disp(['Draw border around: ' str{val}]);
        drawfreehand;
    end
end
