# BNLEARN
Example:
Bài toán: Bệnh nhân bị ảnh hưởng bởi sai khớp cắn loại III (đặc trưng bởi sự nhô ra của cung răng dưới) bị mất cân bằng xương được hình thành từ rất sớm trong đời và điều đó trở nên rõ rệt hơn ở tuổi dậy thì và cho đến khi 
          quá trình trưởng thành của xương hoàn tất. Dự đoán sớm thành công hay thất bại điều trị ở một bệnh nhân Loại III giúp điều chỉnh dễ dàng hơn, nhưng khó thực hiện chỉ từ một số ít các yếu tố quyết định hình thái học 
          là một vấn đề. Lý do cho điều đó là sai khớp cắn Loại III hiếm khi là hậu quả của một bất thường ở một thành phần sọ mặt duy nhất, do đó các phép đo lâm sàng và X quang riêng lẻ có thể ít mang tính biểu thị hơn so 
          với sự tương tác giữa các phép đo.

Tập dữ liệu: bao gồm 143 bệnh nhân với hai bộ số đo ở độ tuổi T1 và T2 (được đo bằng năm) cho các biến sau:
  Treatment: không điều trị (NT), điều trị kết quả xấu (TB), điều trị kết quả tốt (TG).
  Growth: một biến nhị phân có giá trị Tốt hoặc Xấu, được xác định trên cơ sở CoGn-CoA.
  ANB: góc giữa điểm Down A và B (độ).
  IMPA: góc mặt phẳng răng cửa – hàm dưới (độ).
  PPPM: mặt phẳng vòm miệng - góc mặt phẳng hàm dưới (độ).
  CoA: tổng chiều dài hàm trên từ lồi cầu đến điểm Down A (mm).
  GoPg: chiều dài thân hàm dưới từ gonion đến pogonion (mm).
  CoGo: chiều dài nhánh xương hàm dưới từ lồi cầu đến lồi cầu (mm).
  
Công việc:
1.	Tìm hiểu một BN và sử dụng nó để xác định và hình dung sự tương tác giữa các đặc điểm sai khớp cắn Loại III khác nhau trong quá trình phát triển và điều trị.
2.	Kiểm tra tính nhất quán của nó bằng cách xác minh một số giả thuyết được chấp nhận phổ biến về sự tiến triển của sự mất cân bằng xương này.
3.	Cho thấy rằng các đối tượng không được điều trị sẽ phát triển các mô hình tăng trưởng sọ mặt Loại III khác nhau so với các bệnh nhân được điều trị chỉnh nha bằng phương pháp giãn nở hàm trên và điều trị bằng mặt nạ nhanh chóng.
4.	Trong số các bệnh nhân được điều trị, đoạn CoA (chiều dài hàm trên) và góc ANB (mối quan hệ trước-sau của hàm trên và hàm dưới) dường như là các không gian con xương nhận được tác dụng chính của việc điều trị.

file data: prepd-ortho.rda
