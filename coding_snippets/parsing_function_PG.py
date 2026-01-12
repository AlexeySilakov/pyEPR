import wx
import wx.propgrid as pg
import ast
import operator

# -----------------------------
# Safe expression evaluator
# -----------------------------

ALLOWED_OPERATORS = {
    ast.Add: operator.add,
    ast.Sub: operator.sub,
    ast.Mult: operator.mul,
    ast.Div: operator.truediv,
    ast.Pow: operator.pow,
    ast.USub: operator.neg,
}

class SafeEvaluator:
    def __init__(self, variables):
        self.variables = variables

    def eval(self, expr):
        tree = ast.parse(expr, mode='eval')
        return self._eval_node(tree.body)

    def _eval_node(self, node):
        if isinstance(node, ast.Num):
            return node.n

        if isinstance(node, ast.Name):
            if node.id in self.variables:
                return self.variables[node.id]
            raise ValueError(f"Unknown variable: {node.id}")

        if isinstance(node, ast.BinOp):
            op = ALLOWED_OPERATORS[type(node.op)]
            return op(self._eval_node(node.left),
                      self._eval_node(node.right))

        if isinstance(node, ast.UnaryOp):
            op = ALLOWED_OPERATORS[type(node.op)]
            return op(self._eval_node(node.operand))

        raise TypeError("Unsupported expression")


# -----------------------------
# Main Frame
# -----------------------------

class MainFrame(wx.Frame):
    def __init__(self):
        super().__init__(None, title="wxPython PropertyGrid Function Editor",
                         size=(700, 500))

        panel = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)

        # PropertyGrid
        self.pg = pg.PropertyGrid(panel, style=pg.PG_SPLITTER_AUTO_CENTER)
        self.A = {}

        for i in range(1, 11):
            prop = pg.FloatProperty(f"A({i})", value=0.0)
            self.pg.Append(prop)
            self.A[i] = prop

        vbox.Add(self.pg, 1, wx.EXPAND | wx.ALL, 5)

        # Button
        self.btn = wx.Button(panel, label="Define Function")
        self.btn.Bind(wx.EVT_BUTTON, self.on_define_function)
        vbox.Add(self.btn, 0, wx.ALL, 5)

        # Function controls (created later)
        self.func_box = None
        self.spin_B = None

        panel.SetSizer(vbox)

    # -------------------------
    # Button callback
    # -------------------------

    def on_define_function(self, event):
        if self.func_box:
            return  # already created

        parent = self.btn.GetParent()
        sizer = parent.GetSizer()

        hbox = wx.BoxSizer(wx.HORIZONTAL)

        wx.StaticText(parent, label="Function:")
        self.func_box = wx.TextCtrl(parent,
            value="A(1) = A(3) + 0.5 * B",
            size=(300, -1))

        wx.StaticText(parent, label="B:")
        self.spin_B = wx.SpinCtrlDouble(parent, min=-100, max=100,
                                        inc=0.1, initial=0)

        hbox.AddMany([
            (wx.StaticText(parent, label="Function:"), 0, wx.ALL | wx.CENTER, 5),
            (self.func_box, 0, wx.ALL, 5),
            (wx.StaticText(parent, label="B:"), 0, wx.ALL | wx.CENTER, 5),
            (self.spin_B, 0, wx.ALL, 5),
        ])

        sizer.Add(hbox, 0, wx.ALL, 5)
        parent.Layout()

        self.func_box.Bind(wx.EVT_TEXT, self.recompute)
        self.spin_B.Bind(wx.EVT_SPINCTRLDOUBLE, self.recompute)
        self.pg.Bind(pg.EVT_PG_CHANGED, self.recompute)

        self.recompute()

    # -------------------------
    # Core logic
    # -------------------------

    def recompute(self, event=None):
        if not self.func_box:
            return

        expr = self.func_box.GetValue()

        try:
            lhs, rhs = expr.split("=")
            lhs = lhs.strip()
            rhs = rhs.strip()

            if not lhs.startswith("A(") or not lhs.endswith(")"):
                raise ValueError("Left side must be A(n)")

            index = int(lhs[2:-1])

            variables = {
                "B": self.spin_B.GetValue()
            }

            for i in range(1, 11):
                variables[f"A{i}"] = self.pg.GetPropertyValue(f"A({i})")

            # Replace A(n) → An for parser simplicity
            rhs = rhs.replace("A(", "A").replace(")", "")

            evaluator = SafeEvaluator(variables)
            result = evaluator.eval(rhs)

            self.pg.SetPropertyValue(f"A({index})", result)

        except Exception as e:
            # Silent failure by design; easy to add status bar feedback
            pass


# -----------------------------
# App
# -----------------------------

class App(wx.App):
    def OnInit(self):
        frame = MainFrame()
        frame.Show()
        return True


if __name__ == "__main__":
    app = App(False)
    app.MainLoop()
