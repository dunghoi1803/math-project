const { ipcRenderer } = require("electron");
const { sqrt } = require("mathjs");
const math = require("mathjs");

function run() {
let n = document.getElementById('n').value;
let sum_Xi = document.getElementById('sum_Xi').value;
let sum_Zi = document.getElementById('sum_Zi').value;
let sum_Yi = document.getElementById('sum_Yi').value;
let sum_Xi_sqr = document.getElementById('sum_Xi_sqr').value;
let sum_Zi_sqr = document.getElementById('sum_Zi_sqr').value;
let sum_Yi_sqr = document.getElementById('sum_Yi_sqr').value;
let sum_XiZi = document.getElementById('sum_XiZi').value;
let sum_YiZi = document.getElementById('sum_YiZi').value;
let sum_XiYi = document.getElementById('sum_XiYi').value;
let k = document.getElementById('k').value;
let Xo = document.getElementById('Xo').value;
let Zo = document.getElementById('Zo').value;
let hằng_số_t_alpha = document.getElementById('hằng_số_t_alpha').value;

// the first quest

const matrix1 = [[n, sum_Xi ,sum_Zi ], [sum_Xi ,sum_Xi_sqr ,sum_XiZi ], [sum_Zi ,sum_XiZi , sum_Zi_sqr ]]
const matrix_XtX = math.matrix(matrix1)    

const matrix2 = [sum_Yi , sum_XiYi , sum_YiZi ]
const matrix_XtY = math.matrix(matrix2)

const det = math.det(matrix_XtX)

const inverse = (m) => {
    return math.inv(m);
  };
  const matrix_XtX_inverse = inverse(matrix_XtX);

  const betaMu = math.multiply(matrix_XtX_inverse, matrix_XtY)


  // the second quest

let sum_resid_sqr = sum_Yi_sqr - (betaMu.subset(math.index(0))*sum_Yi + betaMu.subset(math.index(1))*sum_XiYi + betaMu.subset(math.index(2))*sum_YiZi);
let sigmaMu_sqr = sum_resid_sqr/(n-k);
let KTC_betaMu = ''
let j = 0
let standardError_betaMu_array = [];

  let variance_betaMu = sigmaMu_sqr*matrix_XtX_inverse.subset(math.index(j, j))
  let standardError_betaMu = Math.sqrt(variance_betaMu)
  standardError_betaMu_array.push(standardError_betaMu)
  let KTC_left = betaMu.subset(math.index(j)) - hằng_số_t_alpha*standardError_betaMu
  let KTC_right = betaMu.subset(math.index(j)) + hằng_số_t_alpha*standardError_betaMu
  KTC_betaMu = '(' + KTC_left + ' ; ' + KTC_right + ')'

  variance_betaMu = sigmaMu_sqr*matrix_XtX_inverse.subset(math.index(j+1, j+1))
  standardError_betaMu = Math.sqrt(variance_betaMu)
  standardError_betaMu_array.push(standardError_betaMu)
  KTC_left = betaMu.subset(math.index(j+1)) - hằng_số_t_alpha*standardError_betaMu
  KTC_right = betaMu.subset(math.index(j+1)) + hằng_số_t_alpha*standardError_betaMu
  KTC_betaMu += ' , (' + KTC_left + ' ; ' + KTC_right + ')'

  variance_betaMu = sigmaMu_sqr*matrix_XtX_inverse.subset(math.index(j+2, j+2))
  standardError_betaMu = Math.sqrt(variance_betaMu)
  standardError_betaMu_array.push(standardError_betaMu)
  KTC_left = betaMu.subset(math.index(j+2)) - hằng_số_t_alpha*standardError_betaMu
  KTC_right = betaMu.subset(math.index(j+2)) + hằng_số_t_alpha*standardError_betaMu
  KTC_betaMu += ' , (' + KTC_left + ' ; ' + KTC_right + ')'


// the third quest

const matrix_Xto = math.matrix([1, Xo, Zo])
const matrix_Xo = math.matrix([[1], [Xo], [Zo]])
let Yo = betaMu.subset(math.index(0)) + betaMu.subset(math.index(1))*Xo + betaMu.subset(math.index(2))*Zo
let variance_YoMu = sigmaMu_sqr*math.multiply(math.multiply(matrix_Xto, matrix_XtX_inverse), matrix_Xo).subset(math.index(0))
let variance_Yo_YoMu = variance_YoMu + sigmaMu_sqr
let standardError_YoMu = Math.sqrt(variance_YoMu)
let standardError_Yo_YoMu = Math.sqrt(variance_Yo_YoMu)
let KTC_left_Yo = Yo - hằng_số_t_alpha*standardError_YoMu;
let KTC_right_Yo = Yo + hằng_số_t_alpha*standardError_YoMu;
let KTC_left_Yo_YoMu = Yo - hằng_số_t_alpha*standardError_Yo_YoMu
let KTC_right_Yo_YoMu = Yo + hằng_số_t_alpha*standardError_Yo_YoMu
const KTC_Y_Xo = KTC_left_Yo+' ; '+ KTC_right_Yo
const KTC_Y_YoMu = KTC_left_Yo_YoMu+' ; '+ KTC_right_Yo_YoMu

// the addtional

let R_sqr = (betaMu.subset(math.index(0))*sum_Yi + betaMu.subset(math.index(1))*sum_XiYi + betaMu.subset(math.index(2))*sum_YiZi - n*((sum_Yi/n)**2))/(sum_Yi_sqr-n*((sum_Yi/n)**2))
let R_sqr_adjusted = 1 - ((1-R_sqr)*(n-1)/(n-k))
let F_tn = (R_sqr/(1-R_sqr))*((n-k)/(k-1))
const T_tn = []
for (i=0;i< standardError_betaMu_array.length;i++) {
  let Ttn = betaMu.subset(math.index(i))/standardError_betaMu_array[i]
  T_tn.push(Ttn)
}

document.getElementById('result1').innerHTML = 'det(XtX)= ' + math.det(matrix_XtX)
document.getElementById('result2').innerHTML = 'Ma trận beta_Mu= ' + betaMu
document.getElementById('result3').innerHTML = 'Tổng bình phương phần dư= ' + sum_resid_sqr
document.getElementById('result4').innerHTML = 'Sigma_Mũ= ' + sigmaMu_sqr**0.5
document.getElementById('result5').innerHTML = 'se(beta_j_Mu)= ' +'[' + standardError_betaMu_array.join(', ') + ']'
document.getElementById('result6').innerHTML = 't_tn= ' +'[' + T_tn.join(', ') + ']'
document.getElementById('result7').innerHTML = 'R^2= ' + R_sqr
document.getElementById('result8').innerHTML = 'R_ngang^2= ' + R_sqr_adjusted
document.getElementById('result9').innerHTML = 'f_tn= ' + F_tn
document.getElementById('result10').innerHTML = 'Khoảng tin cậy beta_j= ' + KTC_betaMu
document.getElementById('result11').innerHTML = 'Y0_Mu= ' + Yo
document.getElementById('result12').innerHTML = 'se(Y0)= ' + standardError_YoMu
document.getElementById('result13').innerHTML = 'se(Y0 - Y0_mu)= ' + standardError_Yo_YoMu 
document.getElementById('result14').innerHTML = 'Dự báo giá trị trung bình của Yo là: ' + '(' + KTC_Y_Xo + ')'
document.getElementById('result15').innerHTML = 'Dự báo giá trị cá biệt của Yo là: ' + '(' + KTC_Y_YoMu + ')'
}
