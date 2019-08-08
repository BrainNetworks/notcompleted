function RenameSheet(FileName,SheetName)

%FileName�����ɵ�excel�����ļ�����
%SheetName�����ɵ�excel����е�sheet�������

% �ж�Excel�Ƿ��Ѿ��򿪣����Ѵ򿪣����ڴ򿪵�Excel�н��в���������ʹ�Excel
try
    % ��Excel�������Ѿ��򿪣���������Excel
    Excel = actxGetRunningServer('Excel.Application');
catch
    % ����һ��Microsoft Excel�����������ؾ��Excel
    Excel = actxserver('Excel.Application');
end

% ����Excel������Ϊ�ɼ�״̬
Excel.Visible = 1;    % set(Excel, 'Visible', 1); 

% �������ļ����ڣ��򿪸ò����ļ��������½�һ���������������棬�ļ���Ϊ����.Excel
if exist(FileName,'file'); 
    Workbook = Excel.Workbooks.Open(FileName);
    % Workbook = invoke(Excel.Workbooks,'Open',filespec_user);
else
    Workbook = Excel.Workbooks.Add;
    % Workbook = invoke(Excel.Workbooks, 'Add');
    Workbook.SaveAs(FileName);
end

% ���ص�ǰ��������
Sheets = Excel.ActiveWorkbook.Sheets;    % Sheets = Workbook.Sheets;
Sheet1 = Sheets.Item(1);    % ���ص�1�������
Sheet1.Activate;    % �����1�����

Sheet1.Name = SheetName;

Workbook.Save   % �����ĵ�
Workbook.Close   % �����ĵ�

% �ر�Excel������
Excel.Visible = 0;    % set(Excel, 'Visible', 0); 
end