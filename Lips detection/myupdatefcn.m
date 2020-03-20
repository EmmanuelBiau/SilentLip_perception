function output_txt = myfunction(~,event_obj)

pos = get(event_obj,'Position');
output_txt = {['X:',num2str(pos(1))],['Y:',num2str(pos(2))]};