function [ramp_length, scaling_factor, Contraction_Type] = dlg_SetTrace(rot_time)

%% initialize text var for the inputs
prompt = {'Enter Rotation duration(s)','Select weeks for scaling the intensity','Contraction Type'};
radioBut_week = {'First part (50%, 75%, 80%, 85%, 90%)', 'Second part (50%, 75%, 80%, 85%, 90%)'};
radioBut_con = {'Concentric', 'Eccentric'};

%% create the object to handle
    % Create a figure window
    fig = figure('Units', 'normalized', 'Position', [0.5 0.5 0.4 0.4], ...
        'Name', 'My Input Dialog', 'NumberTitle', 'off');

    % Create input fields
    uicontrol('Style', 'text', 'String', prompt{1}, ...
        'Units', 'normalized', 'Position', [0.1 0.775 0.3 0.1],'FontSize',14);
    input1Edit = uicontrol('Style', 'edit', 'Units', 'normalized', ...
        'Position', [0.4 0.8 0.3 0.1],'String',rot_time,'Callback', @input1Callback,'FontSize',14);

    % create radio button for weeks intensity
    uicontrol('Style', 'text', 'String', prompt{2}, ...
        'Units', 'normalized', 'Position', [0.1 0.6 0.3 0.1],'FontSize',14);
    radioGroupWk = uibuttongroup('Units', 'normalized', 'Position', [0.4 0.55 0.55 0.2]);
    option1Radio = uicontrol('Style', 'radiobutton', 'String', radioBut_week{1}, ...
        'Units', 'normalized', 'Position', [0.1 0.5 0.8 0.3], 'Parent', radioGroupWk,'FontSize',14);
    option2Radio = uicontrol('Style', 'radiobutton', 'String', radioBut_week{2}, ...
        'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.3], 'Parent', radioGroupWk,'FontSize',14);
    %@input2Callback);

    % Create radio buttons
    uicontrol('Style', 'text', 'String', prompt{3}, ...
        'Units', 'normalized', 'Position', [0.1 0.325 0.3 0.1],'FontSize',14);
    radioGroup = uibuttongroup('Units', 'normalized', 'Position', [0.4 0.3 0.55 0.2]);
    option3Radio = uicontrol('Style', 'radiobutton', 'String', radioBut_con{1}, ...
        'Units', 'normalized', 'Position', [0.1 0.5 0.8 0.3], 'Parent', radioGroup,'FontSize',14);
    option4Radio = uicontrol('Style', 'radiobutton', 'String', radioBut_con{2}, ...
        'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.3], 'Parent', radioGroup,'FontSize',14);

    % Create OK button
    okButton = uicontrol('Style', 'pushbutton', 'String', 'OK', ...
        'Units', 'normalized', 'Position', [0.3 0.05 0.4 0.15], ...
        'Callback', @okButtonCallback);

  %% Initialize output variables
    ramp_length = '';
    scaling_factor = '';
    Contraction_Type = '';

  %% Callback function for Input 1 edit field
    function input1Callback(hObject, ~)
        % Enforce numeric input only
        input = get(hObject, 'String');
        if isempty(input) || isnan(str2double(input))
            set(hObject, 'String', rot_time);
            errordlg('Input must be a number', 'Error', 'modal');
        end
    end

 %% Callback function for Input 2 edit field
    function input2Callback(hObject, ~)
        % Enforce value limits
        input = str2double(get(hObject, 'String'));
        if isnan(input) || input < 0 || input > 100
            set(hObject, 'String', '75');
            errordlg('Input must be a number between 0 and 100', 'Error', 'modal');
        end
    end

  %% Callback function for OK button
    function okButtonCallback(~, ~)
        ramp_length = str2double(get(input1Edit, 'String'));
        %scaling_factor = str2double(get(input2Edit, 'String')) / 100;

        %get training part
        if get(option1Radio, 'Value')
            scaling_factor = [0.5 0.75 0.8 0.85 0.9];
        elseif get(option2Radio, 'Value')
            scaling_factor = [0.5 0.75 0.80 0.85 0.9];
        end


        %get contraction type
        if get(option3Radio, 'Value')
            Contraction_Type = radioBut_con{1};
        elseif get(option4Radio, 'Value')
            Contraction_Type = radioBut_con{2};
        end

        % Close the figure
        close(fig);
    end

    % Wait for the figure to close
    uiwait(fig);

end