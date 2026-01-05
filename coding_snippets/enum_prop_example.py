import wx
import wx.propgrid as wxpg

class PropertyGridFrame(wx.Frame):
    def __init__(self):
        super().__init__(None, title="EnumProperty Choices Example", size=(500, 450))
        
        # Create property grid
        self.pg = wxpg.PropertyGrid(self, style=wxpg.PG_BOLD_MODIFIED)
        
        # Add enum property with initial labels
        initial_labels = ["Red", "Green", "Blue", "Yellow"]
        self.pg.Append(wxpg.EnumProperty("Color", value=1, labels=initial_labels))
        
        # Add another enum property
        status_labels = ["Active", "Inactive", "Pending", "Completed"]
        self.pg.Append(wxpg.EnumProperty("Status", value=2, labels=status_labels))
        
        # Create buttons to demonstrate choice extraction
        button_sizer = wx.BoxSizer(wx.HORIZONTAL)
        
        get_choices_btn = wx.Button(self, label="Get Available Choices")
        get_choices_btn.Bind(wx.EVT_BUTTON, self.on_get_choices)
        button_sizer.Add(get_choices_btn, 0, wx.ALL, 5)
        
        get_selected_btn = wx.Button(self, label="Get Selected Labels")
        get_selected_btn.Bind(wx.EVT_BUTTON, self.on_get_selected)
        button_sizer.Add(get_selected_btn, 0, wx.ALL, 5)
        
        # Add a text control to display results
        self.result_text = wx.TextCtrl(self, style=wx.TE_MULTILINE | wx.TE_READONLY, size=(480, 100))
        
        # Layout
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        main_sizer.Add(self.pg, 1, wx.EXPAND)
        main_sizer.Add(button_sizer, 0, wx.CENTER)
        main_sizer.Add(wx.StaticText(self, label="Results:"), 0, wx.ALL, 5)
        main_sizer.Add(self.result_text, 0, wx.EXPAND | wx.ALL, 5)
        self.SetSizer(main_sizer)
        
    def on_get_choices(self, event):
        """Get the available choices from properties"""
        result = []
        
        # Get Color property
        color_prop = self.pg.GetPropertyByName("Color")
        if color_prop:
            # Get available choices
            choices = color_prop.GetChoices()
            result.append(f"Color choices: {choices}")
            
        # Get Status property
        status_prop = self.pg.GetPropertyByName("Status")
        if status_prop:
            choices = status_prop.GetChoices()
            result.append(f"Status choices: {choices}")
        
        self.result_text.SetValue("\n".join(result))
        
    def on_get_selected(self, event):
        """Get the selected label from properties"""
        result = []
        
        # Get Color property
        color_prop = self.pg.GetPropertyByName("Color")
        if color_prop:
            # Get the selected label
            selected_label = color_prop.GetLabel()
            # Get the integer value
            value = color_prop.GetValue()
            # Get available choices
            choices = color_prop.GetChoices()
            if value < len(choices):
                result.append(f"Color: Selected '{selected_label}' (value: {value}, choice: '{choices[value]}')")
        
        # Get Status property
        status_prop = self.pg.GetPropertyByName("Status")
        if status_prop:
            selected_label = status_prop.GetLabel()
            value = status_prop.GetValue()
            choices = status_prop.GetChoices()
            if value < len(choices):
                result.append(f"Status: Selected '{selected_label}' (value: {value}, choice: '{choices[value]}')")
        
        self.result_text.SetValue("\n".join(result))

app = wx.App()
frame = PropertyGridFrame()
frame.Show()
app.MainLoop()
