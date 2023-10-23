function selection(g, pseudolabel_cats, src, event)

        val = g.Value;
        str = g.String;
        current_selection = str{val};
        disp(['Draw border around: ' str{val}]);
        
        arena_zones(val).(str{val}) = drawpolygon('StripeColor','r');
        assignin('base','arena_zones',arena_zones(val))
%         g.UserData.(current_selection) = drawfreehand;
    end