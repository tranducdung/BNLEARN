library(bnlearn)
library(gplots)
library(MASS)
library(penalized)
library(Rgraphviz)
load("prepd-ortho.rda")


#Tiền xử lý và phân tích dữ liệu thăm dò
#tạo data frame cho tất cả các biến khác nhau và với Growth, Treatment

#Biến Growth và Treatment mang thông tin dư thừa về tiên lượng của bệnh nhân, bằng chứng là sự khác biệt về tỷ lệ bệnh nhân có Tăng trưởng tốt giữa TB và TG.
#Để tránh sự nhầm lẫn có thể xảy ra khi đưa cả hai biến vào mô hình, chúng tôi mã hóa lại Điều trị dưới dạng biến nhị phân trong đó 0 có nghĩa là NT và 1 có nghĩa là TB hoặc TG. Tương tự, chúng ta mã hóa lại Tăng trưởng với 0 nghĩa là Xấu và 1 nghĩa là Tốt.
diff = data.frame(dANB = ortho$ANB2 - ortho$ANB,
                  dPPPM = ortho$PPPM2 - ortho$PPPM,
                  dIMPA = ortho$IMPA2 - ortho$IMPA,
                  dCoA = ortho$CoA2 - ortho$CoA,
                  dGoPg = ortho$GoPg2 - ortho$GoPg,
                  dCoGo = ortho$CoGo2 - ortho$CoGo,
                  dT = ortho$T2 - ortho$T1,
                  Growth = as.numeric(ortho$Growth) - 1,
                  Treatment = as.numeric(ortho$Treatment != "NT"))

#Vì chúng ta sẽ sử dụng BN Gaussian để phân tích nên việc kiểm tra xem các biến có phân phối chuẩn hay không
#, ít nhất là ở mức biên, cũng rất thú vị; và từ các ô bên dưới, điều đó dường như không đúng với tất cả chúng.
par(mfrow = c(2,3), mar = c(4,2,2,2))

for(var in c("dANB", "dPPPM", "dIMPA", "dCoA", "dGoPg", "dCoGo")){
  x = diff[,var]
  hist(x, prob = TRUE, xlab = var, ylab = "", main = "", col = "ivory")
  lines(density(x), lwd = 2, col = "tomato")#Vẽ một đường theo ước tính mật độ các điểm
  curve(dnorm(x, mean = mean(x), sd = sd(x)), from = min(x), to = max(x),
        add = TRUE, lwd = 2, col = "steelblue")#Vẽ đường cong theo hàm mật độ của phân phối chuẩn. 
  
}

#Các biến có được liên kết bởi các mối quan hệ tuyến tính không? Một số trong số đó là như vậy, nhưng không phải tất cả.
pairs(diff[, setdiff(names(diff), c("Growth", "Treatment"))], 
      upper.panel = function(x, y, ...) {
        points(x = x, y = y, col = "grey")
        abline(coef(lm(y ~ x)), col = "tomato", lwd = 2)
        },
      lower.panel = function(x, y, ...) {
      par(usr = c(0, 1, 0, 1))
      text(x = 0.5, y = 0.5, round(cor(x, y), 2), cex = 2)
      }
     ) #so quan hệ tuyến tính theo từng cặp -> chưa biết để làm gì???

#Cuối cùng, chúng ta có thể xem xét liệu các biến có phân cụm theo bất kỳ cách nào hay không vì các biến phân cụm với nhau có nhiều khả năng được liên kết trong BN hơn.
diff.delta = sapply(diff[, 1:6], function(x) x / diff$dT) #lấy giá trị/thời gian
rho = cor(data.frame(diff.delta, Growth = diff$Growth, Treatment = diff$Treatment))#tính hệ số tương quan theo Spearman(tương quan không tham số) Chúng là các hệ số tương quan dựa trên thứ hạng
palette.breaks = seq(0, 1, 0.1)
par(oma = c(2, 2, 2, 1))
heatmap.2(rho, scale = "none", trace = "none", revC = TRUE, breaks = palette.breaks)
#sao lại dùng tính hệ số tương quan theo Spearman

#Chúng ta có thể thấy các cụm trong bản đồ nhiệt: cụm đầu tiên bao gồm dCoGo, dGoPg và dCoA và cụm thứ hai bao gồm Điều trị, dANB và dCoA.
#Cụm đầu tiên rất thú vị về mặt lâm sàng vì nó bao gồm Điều trị và hai biến số đều liên quan đến điểm A của Down, đưa ra một số manh mối về tác dụng chính của việc điều trị.
#-> chưa hiểu sao phân đk 2 cụm như thế này
ug = empty.graph(colnames(rho))
amat(ug) = (rho > 0.4) + 0L - diag(1L, nrow(rho))
graphviz.plot(ug, layout = "fdp", shape = "ellipse")
# tại sao lại lấy > 0.4, chưa hiểu ý nghĩa của đoạn này, có phải lúc nào cũng là -diag(1L, nrow())????

#Mạng bayesian tĩnh
#Learning the Bayesian network
#1.Learning the structure
#b1: làm giảm không gian mô hình tìm được 1 BN tốt hơn. Đưa vào blacklist các cung mã hóa mà t cho là không thể/có thật
#Đưa vào whitelist các cung mà biết là tồn tại.

#action:
#Đưa vào blacklist bất kỳ cung nào trỏ đến dT, Treatment và Growth từ các biến chỉnh nha.
#đưa vào  blacklist vòng cung từ dT đến Treatment. Điều này có nghĩa là việc bệnh nhân có được điều trị hay không không thay đổi theo thời gian.
#đưa vào danh sách đen vòng cung từ Growth đến dT và Treatment. Điều này có nghĩa là việc bệnh nhân có được điều trị hay không không thay đổi theo thời gian và rõ ràng là không thay đổi tùy thuộc vào tiên lượng.

bl = tiers2blacklist(list("dT", "Treatment", "Growth",
                          c("dANB", "dPPPM", "dIMPA", "dCoA", "dGoPg", "dCoGo")))
bl = rbind(bl, c("dT", "Treatment"), c("Treatment", "dT"))
#sao lại nối thêm 2 hàng này???

#Chúng tôi đưa vào danh sách trắng cấu trúc phụ thuộc: dANB → dIMPA ← dPPPM
#Chúng tôi đưa vào danh sách trắng vòng cung từ dT đến Growth để cho phép tiên lượng thay đổi theo thời gian.
wl = matrix(c("dANB", "dIMPA",
              "dPPPM", "dIMPA",
              "dT", "Growth"),
            ncol = 2, byrow = TRUE, dimnames = list(NULL, c("from", "to")))
#chưa hiểu lắm cách chọn này??

#tìm cấu trúc BN tốt hơn bằng cách sử dụng hàm hc() với điểm mặc định (BIC) và toàn bộ khung dữ liệu khác biệt:
dag = hc(diff, whitelist = wl, blacklist = bl)
graphviz.plot(dag, shape = "ellipse", highlight = list(arcs = wl))
#Tuy nhiên, chất lượng của dag chủ yếu phụ thuộc vào việc các biến có được phân phối chuẩn hay không 
#và liệu các mối quan hệ liên kết chúng có tuyến tính hay không; 
#từ phân tích thăm dò, không rõ ràng đó là trường hợp của tất cả chúng. Chúng tôi cũng không biết cung nào thể hiện mối quan hệ bền chặt, nghĩa là chúng có khả năng chống lại sự nhiễu loạn của dữ liệu. Chúng ta có thể giải quyết cả hai vấn đề bằng cách sử dụng boot. Strength() để:
#1. lấy mẫu lại dữ liệu bằng bootstrap;
#2. tìm một mạng riêng biệt từ mỗi mẫu bootstrap;
#3. kiểm tra tần suất mỗi cung có thể xuất hiện trong mạng;
#4. xây dựng một mạng lưới thuận với các vòng cung xuất hiện thường xuyên hơn.
str.diff = boot.strength(diff, R = 200, algorithm = "hc",
                         algorithm.args = list(whitelist = wl, blacklist = bl))
head(str.diff)
#Giá trị trả về của boot. Strength() bao gồm, đối với mỗi cặp nút, cường độ của cung kết nối chúng(Strength) (giả sử tần suất chúng ta quan sát dANB → dPPPM hoặc dPPPM → dANB) 
#và cường độ hướng của nó(direction) (ví dụ: làm thế nào chúng ta thường quan sát dANB → dPPPM khi chúng ta quan sát thấy một cung giữa dANB và dPPPM). 
#boot.Strength() cũng tính toán ngưỡng sẽ được sử dụng để quyết định xem một cung có đủ mạnh để đưa vào mạng đồng thuận hay không.
attr(str.diff, "threshold")# không ra số giống web
#Vì vậy, Averageed.network() lấy tất cả các cung có cường độ ít nhất là attr(str.diff, "threshold") 
#và trả về mạng đồng thuận trung bình, trừ khi chỉ định một ngưỡng khác.
avg.diff = averaged.network(str.diff)
#Vẽ biểu đồ avg.diff bằng Rgraphviz, chúng ta có thể kết hợp thông tin hiện có về độ mạnh của các cung bằng cách sử dụng Strength.plot() thay vì graphviz.plot(). 
#Strength.plot() lấy các đối số tương tự như graphviz.plot() cộng với một ngưỡng và một tập hợp các điểm cắt để xác định cách định dạng từng cung tùy thuộc vào cường độ của nó.
strength.plot(avg.diff, str.diff, shape = "ellipse", highlight = list(arcs = wl))
#So sánh mạng trung bình này với mạng ban đầu
par(mfrow = c(1, 2))
graphviz.compare(avg.diff, dag, shape = "ellipse", main = c("averaged DAG", "single DAG"))
compare(avg.diff, dag, arcs = TRUE)
#Kiểm tra xem cso cung vô hướng không
undirected.arcs(cpdag(dag, wlbl = TRUE))
avg.diff$learning$whitelist = wl
avg.diff$learning$blacklist = bl
undirected.arcs(cpdag(avg.diff, wlbl = TRUE))
#Cuối cùng, chúng ta có thể kết hợp compare() và cpdag() để thực hiện so sánh có nguyên tắc, 
#trong đó chúng ta nói hai cung khác nhau nếu chúng được xác định duy nhất là khác nhau.
compare(cpdag(avg.diff, wlbl = TRUE), cpdag(dag, wlbl = TRUE))
#Một cách xem khác
plot(str.diff)
abline(v = 0.75, col = "tomato", lty = 2, lwd = 2)
abline(v = 0.85, col = "steelblue", lty = 2, lwd = 2)
plot(str.diff)
#tăng thêm ngưỡng và giảm thêm cung(lấy ngưỡng mà ở đó có ít điểm nhất)
nrow(str.diff[str.diff$strength > attr(str.diff, "threshold") &
     str.diff$direction > 0.5, ])
nrow(str.diff[str.diff$strength > 0.75 & str.diff$direction > 0.5, ])
nrow(str.diff[str.diff$strength > 0.85 & str.diff$direction > 0.5, ])
#Mạng đơn giản hơn mà chúng tôi có được bằng cách đặt ngưỡng = 0,85 trong Average.network() được hiển thị bên dưới; 
#chắc chắn là dễ dàng hơn để lập luận từ quan điểm định tính.
avg.simpler = averaged.network(str.diff, threshold = 0.85)
strength.plot(avg.simpler, str.diff, shape = "ellipse", highlight = list(arcs = wl))

#Learning the parameters(Phần này chưa hiểu lắm cần hỏi lại)
#Sau khi học xong cấu trúc, bây giờ chúng ta có thể tìm hiểu các tham số. 
#Vì chúng ta đang làm việc với các biến liên tục nên chúng tôi chọn mô hình hóa chúng bằng GBN. 
#Do đó, nếu chúng ta điều chỉnh các tham số của mạng bằng cách sử dụng maximum likelihood thì chúng ta có rằng mỗi phân phối cục bộ là một hồi quy tuyến tính cổ điển.
fitted.simpler = bn.fit(avg.simpler, diff)
#Chúng ta có thể dễ dàng xác nhận trường hợp đó bằng cách so sánh các mô hình do bn.fit() và lm() tạo ra, chẳng hạn như dANB.
fitted.simpler$dANB
summary(lm(dANB ~ Growth + Treatment, data = diff))
#multivariate Gaussian
mu = rep(0, 3)
R = matrix(c(1, 0.6, 0.5,
             0.6, 1,   0,
             0.5, 0,   1),
             ncol = 3, dimnames = list(c("y", "x1", "x2"), c("y", "x1", "x2")))
#ma trận này tìm theo p2 phân phối chuẩn nhiều chiều? chưa rõ cách tìm cho lắm và sao lại là 3 chiều.
for (rho in seq(from = 0, to = 0.85, by = 0.05)) {
  #cập nhật ma trận tương quan và tạo dữ liệu.
  R[2, 3] = R[3, 2] = rho
  data = as.data.frame(mvrnorm(10000, mu, R)) #Mô phỏng phân phối chuẩn
  #so sánh các mô hình tuyến tính (đầy đủ và tối giản).
  cat("rho:", sprintf("%.2f", rho), "difference in BIC:",
      - 2 * (BIC(lm(y ~ x1 + x2, data = data)) - BIC(lm(y ~ x1, data = data))), "\n")
}
#Nhìn vào BIC này thì đánh giá thế nào nhỉ?

#Nếu ước tính tham số có vấn đề vì bất kỳ lý do gì, chúng ta có thể thay thế chúng bằng một bộ ước tính mới từ một cách tiếp cận khác.
fitted.new = fitted.simpler
fitted.new$dANB = list(coef = c(-1, 2, 2), sd = 1.5)
fitted.new$dANB
fitted.new$dANB = penalized(diff$dANB, penalized = diff[, parents(avg.simpler, "dANB")],
                             lambda2 = 20, model = "linear", trace = FALSE)
fitted.new$dANB

#mẫu chính thức
#Có hai cách tiếp cận chính để xác nhận BN.
#1. Chỉ nhìn vào cấu trúc mạng: nếu mục tiêu chính của việc học BN là xác định các cung và đường dẫn, điều này thường xảy ra khi BN được hiểu là mô hình nhân quả, 
#thì chúng ta có thể thực hiện những gì về cơ bản là phân tích đường dẫn và nghiên cứu cường độ cung.
#2.Nhìn vào BN một cách tổng thể, bao gồm các tham số: nếu mục tiêu chính của việc tìm hiểu BN là sử dụng nó như một mô hình chuyên gia thì chúng ta có thể muốn:
 #- dự đoán giá trị của một hoặc nhiều biến số cho các cá thể mới, dựa trên giá trị của một số biến số khác;
 #- so sánh kết quả của truy vấn CP với kiến thức chuyên môn để xác nhận BN phản ánh kiến thức tốt nhất hiện có về hiện tượng mà chúng ta đang lập mô hình.

#Độ chính xác dự đoán
#Chúng ta có thể đo lường độ chính xác dự đoán của chiến lược học tập đã chọn theo cách thông thường, bằng xác thực chéo. 
#bnlearn cung cấp hàm bn.cv() cho tác vụ này, hàm này thực hiện:
  #- xác thực chéo k-fold
  #- xác thực chéo với các fold do người dùng chỉ định;
  #- giữ lại xác thực chéo
#cho
  #thuật toán học cấu trúc (cấu trúc và các tham số được học từ dữ liệu);
  #thuật toán học tham số (cấu trúc do người dùng cung cấp, tham số được học từ dữ liệu).

#Giá trị trả về của bn.cv() là một đối tượng của lớp bn.kcv (hoặc bn.kcv-list cho nhiều lần chạy xác thực chéo, xem "lớp bn.kcv") có chứa:
  #chỉ mục hàng cho các quan sát được sử dụng làm tập kiểm tra;
  #một đối tượng bn.fit đã học được từ dữ liệu huấn luyện;
  #giá trị của hàm loss;
  #giá trị phù hợp và dự đoán cho loss functions yêu cầu chúng.
#Đầu tiên chúng ta kiểm tra Growth, mã hóa sự tiến triển của sai khớp cắn (0 nghĩa là Xấu và 1 nghĩa là Tốt). Chúng ta kiểm tra chuyển nó trở thành biến rời rạc và tính toán lỗi dự đoán.
xval = bn.cv(diff, bn = "hc", algorithm.args = list(blacklist = bl, whitelist = wl),
             loss = "cor-lw", loss.args = list(target = "Growth", n = 200), runs = 10)
err = numeric(10)

for(i in 1:10) {
  tt = table(unlist(sapply(xval[[i]], '[[', "observed")),
            unlist(sapply(xval[[i]], '[[', "predicted")) > 0.50) # lớn hon 0.5 là lớn hon cái gì
  
  err[i] = (sum(tt) - sum(diag(tt))) / sum(tt)
}

summary(err)

#Các biến khác là liên tục, vì vậy thay vào đó chúng ta có thể ước tính mối tương quan dự đoán của chúng.
predcor = structure(numeric(6),
                    names = c("dCoGo", "dGoPg", "dIMPA", "dCoA", "dPPPM", "dANB"))
for(var in names(predcor)){
  
  xval = bn.cv(diff, bn = "hc", algorithm.args = list(blacklist = bl, whitelist = wl),
              loss = "cor-lw", loss.args = list(target = var, n = 200), runs = 10)
  
  predcor[var] = mean(sapply(xval, function(x) attr(x, "mean")))
  
}
round(predcor, digits = 3)
mean(predcor) # tìm kỳ vọng của các biến

#Khẳng định bằng kiến thức chuyên môn(Đoạn này đọc thì hiểu cái họ check nhưng chưa hiểu đk sao họ lại check những cái này và với các mức xác suất vậy thì mình có kết luận gì)
#Cách khác để xác nhận liệu BN có hợp lý hay không là coi nó như một mô hình hoạt động thật
#và xem liệu nó có thể hiện những phát biểu quan trọng về thế giới mà không được sử dụng làm kiến thức trước đó trong quá trình học hay không. (Nếu không, chúng tôi sẽ chỉ lấy lại thông tin chúng tôi đã đưa vào trước đó!) Một số ví dụ:
  #1.“Sự tăng trưởng quá mức của CoGo sẽ làm giảm PPPM.”
  #Chúng ta kiểm tra giả thuyết này bằng cách tạo các mẫu cho BN được lưu trữ trong fitted.simpler cho cả dCoGo và dPPPM và giả sử không có hoạt động xử lý nào diễn ra. 
  #Khi dCoGo tăng (biểu thị mức tăng trưởng ngày càng nhanh) dPPPM ngày càng trở nên âm (biểu thị sự giảm góc giả sử góc ban đầu là dương.

sim = cpdist(fitted.simpler, nodes = c("dCoGo", "dPPPM"), n = 10^4,
             evidence = (Treatment < 0.5)) #sao lại tạo mẫu ngẫu nhiên từ việc lấy Treatment < 0.5?hya < 0.5 dk cho là không điều trị
plot(sim, col = "grey")
abline(v = 0, col = 2, lty = 2, lwd = 2)
abline(h = 0, col = 2, lty = 2, lwd = 2)
abline(coef(lm(dPPPM ~ dCoGo, data = sim)), lwd = 2)

  #2.Sự tăng trưởng nhỏ của CoGo sẽ tạo ra sự gia tăng trong PPPM
  #Từ hình trên, mức tăng trưởng âm hoặc bằng 0 của CoGo (dCoGo ⋜ 0) tương ứng với mức tăng trưởng dương trong PPPM với xác suất ≈ 0,60. 
   #Thật không may, đối với sự tăng trưởng nhỏ của CoGo (dCoGo ∈ [0, 2]), dPPPM ⋜ 0 với xác suất ≈ 0,50 nên BN không ủng hộ giả thuyết này.
nrow(sim[(sim$dCoGo <= 0) & (sim$dPPPM > 0), ]) / nrow(sim[(sim$dCoGo <= 0), ])
nrow(sim[(sim$dCoGo > 0) & (sim$dCoGo < 2) & (sim$dPPPM > 0), ])
     nrow(sim[(sim$dCoGo) > 0 & (sim$dCoGo < 2),  ])
  #3.Nếu ANB giảm thì IMPA giảm để bù đắp.
  #Kiểm tra bằng mô phỏng như trước đây, chúng ta đang tìm kiếm các giá trị âm của dANB (biểu thị mức giảm giả định góc ban đầu là dương) liên quan đến các giá trị âm của IMPA (tương tự). 
   #Trong hình bên dưới, dANB tỷ lệ thuận với dIMPA, do đó, việc giảm một cái cho thấy cái kia giảm; xu hướng trung bình (đường màu đen) là âm cho cả hai cùng một lúc.
sim = cpdist(fitted.simpler, nodes = c("dIMPA", "dANB"), n = 10^4,
             evidence = (Treatment < 0.5))
plot(sim, col = "grey")
abline(v = 0, col = 2, lty = 2, lwd = 2)
abline(h = 0, col = 2, lty = 2, lwd = 2)
abline(coef(lm(dIMPA ~ dANB, data = sim)), lwd = 2)

  #4.Nếu GoPg tăng mạnh thì cả ANB và IMPA đều giảm. Nếu chúng ta mô phỏng dGoPg, dANB và dIMPA từ BN giả sử dGoPg > 5 
     #(tức là GoPg đang tăng), chúng ta ước tính xác suất dANB > 0 (tức là ANB đang tăng) ở mức ≈ 0,70 và dIMPA < 0 chỉ ở mức ≈ 0,57.
sim = cpdist(fitted.simpler, nodes = c("dGoPg", "dANB", "dIMPA"), n = 10^4,
             evidence = (dGoPg > 5) & (Treatment < 0.5))
nrow(sim[(sim$dGoPg > 5) & (sim$dANB < 0), ]) / nrow(sim[(sim$dGoPg > 5), ])
nrow(sim[(sim$dGoPg > 5) & (sim$dIMPA < 0), ]) / nrow(sim[(sim$dGoPg > 5), ])

  #5.Trị liệu cố gắng ngăn chặn sự giảm ANB. Nếu chúng ta sửa ANB thì có sự khác biệt nào giữa bệnh nhân được điều trị và không được điều trị
  #Đầu tiên, chúng ta có thể kiểm tra mối quan hệ giữa Treatment và growth đối với những bệnh nhân có dANB ≈ 0 mà không cần bất kỳ sự can thiệp nào (sử dụng BN mà chúng tôi đã học được từ dữ liệu).
sim = cpdist(fitted.simpler, nodes = c("Treatment", "Growth"), n = 5 * 10^4,
              evidence = abs(dANB) < 0.1)
tab = table(TREATMENT = sim$Treatment < 0.5, GOOD.GROWTH = sim$Growth > 0.5)
round(prop.table(tab, margin = 1), 2)
  #P(GOOD.GROWTH ∣ TREATMENT) ước tính là khác nhau đối với bệnh nhân được điều trị và không được điều trị (≈ 0,65 so với ≈ 0,52).
  #Nếu chúng ta mô phỏng một can thiệp chính thức (a la Judea Pearl) và đặt dANB = 0 bên ngoài (do đó làm cho nó độc lập với parents của nó và loại bỏ các cung tương ứng)
  #chúng ta có rằng GOOD.GROWTH thực tế có phân phối giống nhau cho cả bệnh nhân được điều trị và không được điều trị và do đó trở nên độc lập với TREATMENT
  #Điều này cho thấy tiên lượng thuận lợi thực sự được xác định bằng cách ngăn chặn những thay đổi trong ANB và các thành phần khác của điều trị (nếu có) sau đó sẽ trở nên không liên quan.
avg.mutilated = mutilated(avg.simpler, evidence = list(dANB = 0))
fitted.mutilated = bn.fit(avg.mutilated, diff)
fitted.mutilated$dANB = list(coef = c("(Intercept)" = 0), sd = 0)
sim = cpdist(fitted.mutilated, nodes = c("Treatment", "Growth"), n = 5 * 10^4,
             evidence = TRUE)
tab = table(TREATMENT = sim$Treatment < 0.5, GOOD.GROWTH = sim$Growth > 0.5)
round(prop.table(tab, margin = 1), 2)
  #6. Trị liệu cố gắng ngăn chặn sự giảm ANB. Nếu chúng ta khắc phục ANB thì có sự khác biệt nào giữa bệnh nhân được điều trị và không được điều trị không?
  #Một cách để đánh giá điều này là kiểm tra xem góc giữa điểm A và điểm B (ANB) có thay đổi giữa bệnh nhân được điều trị và không được điều trị trong khi vẫn giữ GoPg cố định hay không.
sim.GoPg = cpdist(fitted.simpler, nodes = c("Treatment", "dANB", "dGoPg"),
                  evidence = abs(dGoPg) < 0.1)
  #Giả sử GoPg không thay đổi, góc giữa điểm A và điểm B tăng đối với bệnh nhân được điều trị (giá trị âm mạnh biểu thị sự mất cân bằng theo chiều ngang, do đó tỷ lệ thay đổi dương cho thấy mức độ mất cân bằng giảm) và giảm đối với bệnh nhân không được điều trị (mất cân bằng dần dần xấu đi theo thời gian).
sim.GoPg$Treatment = c("UNTREATED", "TREATED")[(sim.GoPg$Treatment > 0.5) + 1L]
mean(sim.GoPg[sim.GoPg$Treatment == "UNTREATED", "dANB"])
mean(sim.GoPg[sim.GoPg$Treatment == "TREATED", "dANB"])
boxplot(dANB ~ Treatment, data = sim.GoPg)

#2.Mô hình số 2: mạng Bayes động:
#là các BN mô hình hóa các quá trình ngẫu nhiên: mỗi biến được liên kết với một nút khác nhau trong mỗi thời điểm được mô hình hóa. 
#(Thông thường, chúng tôi giả định rằng quy trình là Markov cấp một, vì vậy chúng tôi có hai điểm thời gian trong BN: t và t - 1.) Tuy nhiên, chúng tôi khám phá nó nhằm mục đích minh họa cách học và sử dụng một BN như vậy trong bnlearn.
#Dữ liệu chúng ta sử dụng cho mô hình này là dữ liệu mà chúng ta lưu trữ vào ortho khi bắt đầu phân tích. 
#Tuy nhiên, chúng ta chọn sử dụngTreatment thay vì GROWTH làm biến số để thể hiện thực tế rằng các đối tượng có thể đang tiến hành điều trị y tế. 
#Nguyên nhân là vì GROWTH là một biến đo lường tiên lượng tại thời điểm đo lần thứ hai và không xác định được giá trị của nó tại thời điểm đo lần đầu; trong khi Treatment ở cả hai thời điểm đều giống nhau.

#Học cấu trúc
#Đầu tiên, chúng ta chia các biến thành ba nhóm: các biến tại thời điểm t2, các biến tại thời điểm t1 = t2 - 1 và các biến không phụ thuộc vào thời gian vì chúng có cùng giá trị tại t1
const = "Treatment"
t2.variables = grep("2$", names(ortho), value = TRUE)
t1.variables = setdiff(names(ortho), c(t2.variables, const))
#Sau đó, chúng tôi có một blacklist trong đó:
 #- đưa vào blacklist tất cả các cung từ các biến lâm sàng đến T1, T2 và Treatment vì chúng ta biết rằng độ tuổi và chế độ điều trị không bị quyết định bởi các phép đo lâm sàng.
 #- đưa vào blacklist tất cả các cung đi vào Treatment và vào tất cả các biến tại thời điểm t1, bởi vì chúng tôi giả định rằng các cung giữa các biến tại thời điểm t1 giống với các biến tương ứng ở thời điểm t2 và việc tìm hiểu chúng hai lần là vô nghĩa.
 #- đưa vào blacklist tất cả các cung từ t2 đến t1.
roots = expand.grid(from = setdiff(names(ortho), c("T1", "T2", "Treatment")),
                    to = c("T1", "T2", "Treatment"), stringsAsFactors = FALSE)
empty.t1 = expand.grid(from = c(const, t1.variables), to = c(const, t1.variables),
                       stringsAsFactors = FALSE)
bl = rbind(tiers2blacklist(list(t1.variables, t2.variables)), roots, empty.t1)
#Ngược lại, chỉ đưa vào danh sách trắng cung T1 → T2, vì độ tuổi ở lần đo thứ hai rõ ràng phụ thuộc vào độ tuổi ở lần đo đầu tiên.
wl = data.frame(from = c("T1"), to = c("T2"))
#Cuối cùng chúng ta có thể tìm hiểu cấu trúc của BN với bl và wl.
dyn.dag = tabu(ortho, blacklist = bl, whitelist = wl)
dyn.dag

gR = graphviz.plot(dyn.dag, shape = "rectangle", render = FALSE)
sg0 = list(graph = subGraph(const, gR), cluster = TRUE)
sg1 = list(graph = subGraph(t1.variables, gR), cluster = TRUE)
sg2 = list(graph = subGraph(t2.variables, gR), cluster = TRUE)
gR = layoutGraph(gR, subGList = list(sg0, sg1, sg2),
                 attrs = list(graph = list(rankdir = "LR")))
nodeRenderInfo(gR)$fill[t1.variables] = "tomato"
nodeRenderInfo(gR)$fill[t2.variables] = "gold"
renderGraph(gR)
#Như trong mô hình trước, phương pháp điều trị tác động lên ANB: các cung duy nhất đi ra khỏi Điều trị là Điều trị → ANB2 và Điều trị → CoA2. Một lần nữa cả hai nút con đều liên quan đến điểm A của Down.

#Model averaging in structure learning
#Chúng ta muốn đánh giá tính ổn định của cấu trúc BN động này giống như chúng tôi đã làm với BN tĩnh trước đó và chúng ta có thể thực hiện lại điều đó với boot. Strength() và averated.network().
dyn.str = boot.strength(ortho, R = 200, algorithm = "tabu",
                        algorithm.args = list(blacklist = bl, whitelist = wl))
plot(dyn.str)

dyn.avg = averaged.network(dyn.str)
dyn.avg
#tính trung bình Dyn.avg và dyn.dag gần như giống hệt nhau: chúng chỉ khác nhau bởi hai cung. Điều này cho thấy rằng việc học cấu trúc tạo ra đầu ra ổn định.
#unlist(compare(dyn.dag, dyn.avg)) ko chạy đk

par(mfrow = c(1, 2))
graphviz.compare(dyn.dag, dyn.avg, shape = "rectangle")

#Learning the parameters
#Vì Treatment là một biến rời rạc nên BN là CLGBN. Điều này có nghĩa là các nút liên tục có Phương pháp xử lý là nút gốc có thông số hóa khác với các nút còn lại.
dyn.fitted = bn.fit(dyn.avg, data = ortho)
dyn.fitted$ANB2
#Như chúng ta có thể thấy, ANB2 phụ thuộc vào ANB (vì cùng một biến ở thời điểm trước đó) và Treatment. ANB là liên tục nên nó được sử dụng làm biến hồi quy cho ANB2. Việc xử lý là rời rạc và xác định các thành phần của hỗn hợp hồi quy tuyến tính.

#Model validation and inference
 #1.ANB thay đổi như thế nào từ lần đo thứ nhất sang lần đo thứ hai với các chế độ điều trị khác nhau?
   #Chúng ta có thể tạo các cặp (ANB, ANB2) với cpdist() có điều kiện Điều trị bằng NT, TB và TG và xem xét sự phân bổ của chúng.
nt = cpdist(dyn.fitted, nodes = c("ANB", "ANB2"), evidence = (Treatment == "NT"))
tb = cpdist(dyn.fitted, nodes = c("ANB", "ANB2"), evidence = (Treatment == "TB"))
tg = cpdist(dyn.fitted, nodes = c("ANB", "ANB2"), evidence = (Treatment == "TG"))
effect = data.frame(diff = c(nt[, 2] - nt[, 1], tb[, 2] - tb[, 1], tg[, 2] - tg[, 1]),
                    treatment = c(rep("NT", nrow(nt)), rep("TB", nrow(tb)), rep("TG", nrow(tg)))
)
by(effect$diff, effect$treatment, FUN = mean) #tính trung bình diff, index treatment

col = c("steelblue", "gold", "tomato")
lattice::densityplot(~ diff, groups = treatment, data = effect, col = col, lwd = 2, bw = 2, ylim = c(0, 0.20),
                     key = list(text = list(c("untreated", "treated with bad results", "treated with good results")),
                                col = col, lines = TRUE, corner = c(0.98, 0.98), lwd = 2))
#Chúng ta biết rằng liệu pháp cố gắng ngăn chặn sự giảm ANB; và điều này phù hợp với thực tế là phân phối của NT nằm ở bên trái của TB nằm ở bên trái của TG. 
#Tình trạng bệnh nhân không được điều trị tiếp tục xấu đi; những bệnh nhân được điều trị không hiệu quả không thực sự cải thiện nhưng tình trạng của họ cũng không xấu đi; và bệnh nhân được điều trị có hiệu quả cải thiện.

  #2. Sự phát triển của ANB trông như thế nào đối với các chế độ điều trị khác nhau khi bệnh nhân già đi?
  #Giả sử điều kiện ban đầu của ANB bằng 1 ở 5 tuổi, chúng ta có thể dự đoán lặp đi lặp lại ANB2 cho độ tuổi hiện tại + 3 tuổi để xây dựng quỹ đạo từ thời thơ ấu đến tuổi trưởng thành.
  #Nhưng điều này nêu bật một trong những hạn chế chính của mô hình này: giả định rằng các phụ thuộc xác suất là tuyến tính có nghĩa là quỹ đạo của ANB2 cũng sẽ gần như tuyến tính. 
  #Điều đó là không thực tế: chúng ta sẽ ngừng điều trị trước khi tạo ra sự mất cân bằng theo hướng khác và quá trình tăng trưởng chắc chắn sẽ tác động đến sự phát triển của xương theo cách phi tuyến tính.
#Giả định được điều trị
intervals = data.frame(T1   = c(5, 8, 11, 14, 17),
                       T2   = c(8, 11, 14, 17, 20),
                       ANB  = c(-1, NA, NA, NA, NA),
                       ANB2 = c(NA, NA, NA, NA, NA)
)

for (i in seq(nrow(intervals))){
  
  predictor = data.frame(Treatment = factor("TG", levels = c("NT", "TB", "TG")),
                         T1 = intervals[i, "T1"],
                         T2 = intervals[i, "T2"],
                         ANB = intervals[i, "ANB"]
  )
  intervals[i, "ANB2"] = predict(dyn.fitted, node = "ANB2", data = predictor,
                                 method = "bayes-lw", from = names(predictor), n = 1000)
  if (i < nrow(intervals))
    intervals[i + 1, "ANB"] = intervals[i, "ANB2"]

}
print(intervals)
#Giả định không được điều trị
intervals2 = data.frame(T1  = c(5, 8, 11, 14, 17),
                       T2   = c(8, 11, 14, 17, 20),
                       ANB  = c(-1, NA, NA, NA, NA),
                       ANB2 = c(NA, NA, NA, NA, NA)
)
for (i in seq(nrow(intervals2))){
  
  predictor = data.frame(Treatment = factor("NT", levels = c("NT", "TB", "TG")),
                         T1 = intervals2[i, "T1"],
                         T2 = intervals2[i, "T2"],
                         ANB = intervals2[i, "ANB"]
  )
  intervals2[i, "ANB2"] = predict(dyn.fitted, node = "ANB2", data = predictor,
                                 method = "bayes-lw", from = names(predictor), n = 1000)
  if (i < nrow(intervals2))
    intervals2[i + 1, "ANB"] = intervals2[i, "ANB2"]
  
}
print(intervals2)

#Quỹ đạo mô phỏng của CoA thực tế hơn: nó chậm lại theo tuổi tác. Điều này không giống ANB và nó xảy ra vì CoA2 phụ thuộc vào cả T1 và T2. (ANB2 không phụ thuộc vào cái nào cả.)
intervals = data.frame(T1   = c(5, 8, 11, 14, 17),
                       T2   = c(8, 11, 14, 17, 20),
                       ANB  = c(-1, NA, NA, NA, NA),
                       ANB2 = c(NA, NA, NA, NA, NA),
                       CoA  = c(75, NA, NA, NA, NA),
                       CoA2 = c(NA, NA, NA, NA, NA)
)

for (i in seq(nrow(intervals))) {
  
  predictor = data.frame(Treatment = factor("TG", levels = c("NT", "TB", "TG")),
                         T1 = intervals[i, "T1"],
                         T2 = intervals[i, "T2"],
                         ANB = intervals[i, "ANB"],
                         CoA = intervals[i, "CoA"]
  )
  
  #để thực hiện dự đoán chung, hiện tại không thể thực hiện được với predict().
  dist = cpdist(dyn.fitted, nodes = c("ANB2", "CoA2"),
                evidence = as.list(predictor), method = "lw")

  weights = attr(dist, "weights")
  weights
  intervals[i, "ANB2"] = weighted.mean(dist$ANB2, weights)
  intervals[i, "CoA2"] = weighted.mean(dist$CoA2, weights)
  
  if (i < nrow(intervals)) {
    
    intervals[i + 1, "ANB"] = intervals[i, "ANB2"]
    intervals[i + 1, "CoA"] = intervals[i, "CoA2"]
    
   }
}

